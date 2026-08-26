# 09 — Maintaining it

> **Concept** → The code → Do it yourself → Check yourself

Most of a project's life is maintenance. This chapter is about the parts that
decide whether that life is pleasant.

## Continuous integration

`.github/workflows/history-timeline-ci.yml`. On every push: typecheck → migrate →
seed → test → build.

```yaml
services:
  postgres:
    image: postgres:16-alpine
env:
  TEST_DATABASE_URL: postgres://historian:historian@localhost:5432/history_timeline
```

Because `TEST_DATABASE_URL` is set, **CI runs the Postgres contract tests that
are skipped locally.** Fast loop on your laptop, full confidence in CI — nothing
is permanently skipped.

Four details worth stealing:

- **`npm ci`, not `npm install`.** `ci` installs exactly the lockfile. `install`
  may resolve a newer patch version, so CI tests something you never ran.
- **`concurrency: cancel-in-progress`.** A new push cancels the previous run. No
  point verifying a commit that has been replaced.
- **`paths:` filter.** The workflow only runs when this project changes.
- **Build after test.** A green suite that does not compile helps nobody.

**Make CI required to merge.** A check people can ignore is decoration.

## Logging

Rules that survive contact with a real incident:

**Log at the boundaries, not everywhere.** Requests in, errors out. Logging every
function entry produces noise that hides the one line that mattered.

**Never log secrets.** Passwords, tokens, connection strings, personal data. Once
it is in CloudWatch it is in your compliance scope.

**Log the unexpected in full, tell the user nothing:**

```ts
console.error('[api] unhandled error', error);          // everything, for us
res.status(500).json({ error: { code: 'INTERNAL_ERROR',
  message: isProduction ? 'Something went wrong.' : ... } });   // nothing, for them
```

**Structured logs when you graduate.** `console.log` is fine here. In production
JSON logs are queryable:

```ts
console.log(JSON.stringify({ level: 'info', event: 'request',
  method: req.method, path: req.path, status: res.statusCode, durationMs, requestId }));
```

CloudWatch Logs Insights can then answer "p99 latency on `/api/timeline` last
Tuesday" — which no amount of `console.log('here')` will.

A **request id** on every log line, propagated through the request, is what turns
"an error happened" into "here is that user's whole journey".

## Dependencies

```bash
npm outdated              # what has moved
npm audit                 # known vulnerabilities
npm audit fix             # the safe subset
```

**Update regularly and in small batches.** Twelve months of deferred updates is
not one big job; it is a wall of simultaneous breaking changes with no way to
tell which one broke you.

Semver, and how much to trust it:

- `~4.18.2` — patch only. Bug fixes.
- `^4.18.2` — minor + patch. New features, no breaking changes *by promise*.
- `4.18.2` — exactly this.

The lockfile (`package-lock.json`) is what actually makes builds reproducible.
**Commit it.**

Before adding a dependency, ask: how many lines would this take me to write? Is
it maintained? How many transitive dependencies does it drag in? This project
uses five runtime dependencies on the server and two on the client — the
migrator was ~30 lines rather than a library, and now you understand migrations.

## Configuration

Validated once, at boot, in `Config.fromEnv()`:

```ts
const parsed = EnvSchema.safeParse(env);
if (!parsed.success) throw new Error(`Invalid environment configuration:\n${issues}`);
```

**A misconfigured server should fail loudly on startup, not mysteriously on the
first request at 3am.**

`.env.example` is committed; `.env` is git-ignored. In production, config comes
from AWS Secrets Manager (chapter 10) — never from a file in the image.

## Where a bug gets fixed

The most useful maintenance habit in this codebase. When something is wrong, ask
**which layer owns this rule?**

| Symptom | Fix belongs in |
|---|---|
| Year 3000 is accepted | `domain/Year.ts` |
| Empty country list where events exist | `services/TimelineService.ts` |
| Wrong HTTP status | `api/errorHandler.ts` or the controller |
| Wrong sort order | the repositories — **and the contract test** |
| Events shown for the wrong year | `App.tsx` state, or `useAsyncResource` |
| An uncited event got in | `domain/HistoricalEvent.ts` — and ask how it got past |

Fixing a symptom in the wrong layer is how a codebase rots: the rule ends up
duplicated in three places, and the next person changes two of them.

## Code review

What to actually look for, roughly in order:

1. **Does it do what it says?** Read the tests first — they state the intent.
2. **Is the rule in the right layer?** (The table above.)
3. **What happens on the unhappy path?** Empty, missing, malformed, slow.
4. **Is it testable?** If it is hard to test, it is usually too coupled.
5. **Will the next person understand it?** Comments explain *why*; the code says
   *what*.

Not worth a comment: formatting (a formatter's job) or style preferences with no
behavioural difference. Reviews that bikeshed on naming while missing an N+1
query are worse than no review.

## Documentation that stays true

Documentation lies as soon as the code changes. The defence is to **make as much
of it executable as possible**:

- Tests document behaviour and fail when it changes.
- Types document shape and fail when it changes.
- Comments explain *why* — the part that does not go stale, because reasons
  outlive implementations.

Compare:

```ts
// ❌ says what the next line says, and will outlive its truth
// Set the year
this.year = year;

// ✅ explains a decision the code cannot
// COUNT(*) returns BIGINT, which `pg` hands back as a STRING to avoid silently
// truncating values above 2^53. Parse it explicitly.
```

## Performance, when it is time

**Measure first.** The two things most likely to matter here:

- **N+1 queries** — a query inside a loop. `countByYearRange` exists to avoid
  one.
- **Missing indexes** — `EXPLAIN ANALYZE` on your slowest query.

Then caching: `/api/timeline` changes rarely and would happily sit behind a
`Cache-Control: public, max-age=300` header and CloudFront.

Do not optimise before measuring. The bottleneck is rarely where you think.

## Do it yourself

1. **Read the CI file end to end.** Explain what each step protects against.
2. **Make CI fail.** Push a branch with a deliberate type error and watch which
   step catches it.
3. **Add request logging** middleware: method, path, status, duration, in JSON.
4. **Audit the dependencies.** `npm outdated`, then update one minor version and
   run the suite.

## Check yourself

- Why does CI run the Postgres contract tests when local `npm test` skips them?
- Why validate configuration at startup rather than on first use?
- What makes a comment worth writing?

<details>
<summary>Answers</summary>

- **CI runs them**: the local loop optimises for speed (no database required);
  CI optimises for confidence and has a Postgres service container. Setting
  `TEST_DATABASE_URL` in CI means the skip is a local convenience, not a
  permanent hole.
- **Startup validation**: a config error at boot is one obvious failure before
  any traffic arrives. The same error on first use is an outage at 3am, probably
  reported as "the site is broken sometimes".
- **Worth writing**: it explains *why* — a decision, a trade-off, a non-obvious
  constraint, a bug that motivated the shape of the code. Comments that restate
  the next line are noise that eventually becomes a lie.

</details>

→ Next: [10 — Deploying to AWS](10-aws-deployment.md)
