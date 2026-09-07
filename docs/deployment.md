# Deploying the read-only service

Only **one** component of G4-WATCH may face a network.

| Component | Binds | Authentication | Executes stages | Publishable |
|---|---|---|---|---|
| `web.runner` (operator console) | `127.0.0.1:8010` | **none** | yes | **never** |
| `web.backend` (read-only service) | `0.0.0.0:8000` | configurable | no | yes, behind the proxy |

The console launches pipeline stages as subprocesses and has no
authentication of any kind. Exposing it to a network is remote code
execution on the analysis host. It has no deployment story on purpose:
run it over SSH on the machine that holds the corpus.

## 1. Choose an authentication backend

The default is `OpenAccessBackend` — no login. That is honest for a
service publishing already-public pipeline artifacts on a trusted
network, and wrong for anything else.

```python
from web.backend.app import create_app
from web.backend.auth import LocalUserBackend, hash_password

users = {
    "asharma": (hash_password("..."), "A. Sharma", frozenset({"viewer"})),
}
app = create_app(LocalUserBackend(users))
```

`hash_password` uses PBKDF2-HMAC-SHA256 from the standard library, and the
hash carries its own salt and iteration count so the count can be raised
later without invalidating existing credentials. Never store a plain
password, and never a bare SHA-256 of one.

For institutional SSO, implement `AuthBackend` — that is what the abstract
class is for. Nothing else needs to change.

## 2. Set the session secret

A backend that requires a login will **refuse to start** without one:

```bash
export G4WATCH_SESSION_SECRET="$(python -c 'import secrets; print(secrets.token_urlsafe(48))')"
```

At least 32 bytes. The refusal is deliberate — a secret generated per
process would log every user out on restart and, with more than one
worker behind the proxy, would reject its own cookies at random.

Keep it out of the repository and out of the image. A file read by the
unit, or the orchestrator's secret store.

## 3. Terminate TLS at the proxy

`web/proxy/nginx.conf` expects:

```
/etc/nginx/certs/fullchain.pem
/etc/nginx/certs/privkey.pem
```

It will not start without them. That is intentional: a proxy that quietly
fell back to plain HTTP would transmit the session cookie in clear, and
the cookie's `Secure` flag would silently stop it working rather than
warn you.

The config sets HSTS for one year **without** `preload`. Raise it to
preload only once certificate renewal is proven — a preload entry is
effectively irreversible for the domain.

Rate limits: `/login` at 5 requests/minute per IP, everything else at
30/s. The login limit is the one that matters, because there request
throughput is guess throughput.

### Running the app on plain HTTP for development

```bash
export G4WATCH_INSECURE_COOKIES=1
```

Opt-out, never the default. Do not set it in a deployment.

## 4. Keep the audit log

Every command the console runs is appended to
`results/audit/console_audit.tsv`: timestamp, event, command, pathogen,
the exact argv, exit code, duration, user.

Put it somewhere durable:

```bash
export G4WATCH_AUDIT_LOG=/var/lib/g4watch/console_audit.tsv
```

This is the operational record. It is **not** the scientific record —
that is `data/atlases/testing_ledger.tsv`, which the D.H1 gate reads.
Nothing in the audit log can authorise a scoring run, and the gate never
consults it.

## 5. Back up the testing ledger

Everything else in `data/` is regenerable from
`(git commit, pathogen config, accession list)`. The ledger is not: it
records what was tested and what the verdict was, and the gate reads it
rather than trusting a configuration flag.

Back it up independently of the regenerable artifacts, and test a restore.
A restore must never rewrite history — the file is append-only.

## What is still missing

- No images have been pushed to a registry, so `nextflow.config` pins by
  tag rather than digest.
- No Prometheus metrics endpoint.
- No alerting on run failure.
