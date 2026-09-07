"""The reverse proxy is the TLS boundary for the only publishable service.

Static assertions over nginx.conf. They do not need nginx: what matters
is that the directives are present and correct, and a missing HSTS header
or an absent login rate limit is a silent weakening that no functional
test would notice.
"""

from __future__ import annotations

import re
from pathlib import Path

CONF = Path("web/proxy/nginx.conf").read_text()

#: Directives only. A config that explains in a comment why it does NOT
#: enable something would otherwise be flagged by a search for that word —
#: the fix and the defect look identical to grep.
DIRECTIVES = "\n".join(
    line for line in CONF.splitlines() if not line.lstrip().startswith("#")
)


def test_plain_http_redirects_rather_than_serving():
    """A session cookie is marked Secure, so serving over HTTP would
    silently break login instead of protecting it."""
    assert re.search(r"return 30[18] https://", CONF), "port 80 does not redirect to HTTPS"


def test_acme_renewal_stays_reachable_over_http():
    """A certificate that cannot renew is a scheduled outage."""
    assert "/.well-known/acme-challenge/" in CONF


def test_tls_is_configured():
    assert "listen 443 ssl" in CONF
    assert "ssl_certificate" in CONF and "ssl_certificate_key" in CONF


def test_only_modern_tls_versions_are_offered():
    protocols = re.search(r"ssl_protocols\s+([^;]+);", CONF).group(1).split()
    assert set(protocols) <= {"TLSv1.2", "TLSv1.3"}, f"weak protocol offered: {protocols}"


def test_hsts_is_set_with_a_meaningful_lifetime():
    match = re.search(r"Strict-Transport-Security\s+\"max-age=(\d+)", CONF)
    assert match, "no HSTS header"
    assert int(match.group(1)) >= 15_552_000, "HSTS lifetime under six months"


def test_hsts_does_not_preload_before_renewal_is_proven():
    """A preload entry is effectively irreversible for the domain."""
    assert "preload" not in DIRECTIVES


def test_the_login_endpoint_is_rate_limited_far_harder_than_the_rest():
    """It is the one endpoint where request throughput is guess throughput."""
    zones = dict(re.findall(r"limit_req_zone[^;]*zone=(\w+):[^;]*rate=(\S+?);", CONF))
    assert "g4_login" in zones, "no login rate-limit zone"
    login_rate = zones["g4_login"]
    assert login_rate.endswith("r/m"), f"login limited per second, not per minute: {login_rate}"
    assert int(login_rate.removesuffix("r/m")) <= 10


def test_the_login_location_actually_applies_that_zone():
    """Declaring a zone and never using it is the easy mistake."""
    block = re.search(r"location = /login \{(.*?)\n    \}", CONF, re.S)
    assert block, "no dedicated /login location"
    assert "limit_req zone=g4_login" in block.group(1)


def test_security_headers_are_present():
    for header in (
        "X-Content-Type-Options",
        "X-Frame-Options",
        "Referrer-Policy",
        "Content-Security-Policy",
    ):
        assert header in CONF, f"{header} is not set"


def test_the_content_security_policy_forbids_inline_script():
    policy = re.search(r"Content-Security-Policy \"([^\"]+)\"", CONF).group(1)
    script = re.search(r"script-src ([^;]+)", policy).group(1)
    assert "'unsafe-inline'" not in script and "'unsafe-eval'" not in script


def test_the_server_version_is_not_advertised():
    assert "server_tokens off" in CONF


def test_upload_size_is_capped_on_a_read_only_service():
    assert re.search(r"client_max_body_size\s+\d+m", CONF)
