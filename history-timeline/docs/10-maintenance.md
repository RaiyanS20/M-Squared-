# 10 — Maintaining it

> **Concept** → The code → Do it yourself → Check yourself

Most of a project's life is maintenance. This chapter is about the parts that
decide whether that life is pleasant.

## Continuous integration

`.github/workflows/history-timeline-ci.yml`. On every push: install → migrate →
seed → test → boot the server and prove it serves.

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
- **`concurrency: cancel-in-progress`.** A new push cancels the previous run.
- **`paths:` filter.** The workflow only runs when this project changes.
- **A smoke test instead of a build step.** There is nothing to compile here, so
  CI starts the server and curls it instead — proving the thing actually boots
  and serves both the API and the page.

**Make CI required to merge.** A check people can ignore is decoration.

## What replaces a compiler

This stack has no type checker, so a few habits carry that weight:

1. **Validate at every boundary** — the domain constructors and `Config.fromEnv`.
2. **Contract tests** for anything with more than one implementation.
3. **A linter.** This project ships without one to keep the dependency list
   honest, but adding ESLint is the single highest-value tool you can add:
   ```bash
   npm install --save-dev eslint
   npx eslint --init
   ```
   It catches unused variables, unreachable code, `==` vs `===`, and shadowed
   names — a real slice of what a compiler would.
4. **When the project grows past what tests comfortably cover, add TypeScript.**
   The architecture is ready for it: the layers, the constructors and the
   contracts are all where a type checker would want them. Adopting it is
   incremental — you can start with `// @ts-check` and JSDoc comments in single
   files, with no build step at all.

## Logging

Rules that survive contact with a real incident:

**Log at the boundaries, not everywhere.** Requests in, errors out. Logging every
function entry produces noise that hides the one line that mattered.

**Never log secrets.** Passwords, tokens, connection strings, personal data.

**Log the unexpected in full, tell the user nothing:**

```js
console.error('[api] unhandled error', error);              // everything, for us
res.status(500).json({ error: { code: 'INTERNAL_ERROR',
  message: isProduction ? 'Something went wrong.' : ... } });  // nothing, for them
```

**Structured logs when you graduate.** `console.log` is fine here. In production,
JSON logs are queryable:

```js
console.log(JSON.stringify({ level: 'info', event: 'request',
  method: req.method, path: req.path, status: res.statusCode, durationMs, requestId }));
```

CloudWatch Logs Insights can then answer "p99 latency on `/api/timeline` last
Tuesday", which no amount of `console.log('here')` will. A **request id** on
every line turns "an error happened" into "here is that user's whole journey".

## Dependencies

This project has **four runtime dependencies** (`express`, `cors`, `pg`,
`dotenv`) and **one dev dependency** (`jsdom`). That is deliberate, and it is
worth keeping.

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
it maintained? How many transitive dependencies does it drag in? The migrator
here is ~30 lines rather than a library, and now you understand migrations.

## Configuration

Validated once, at boot, in `Config.fromEnv()`. **A misconfigured server should
fail loudly on startup, not mysteriously on the first request at 3am.**

`.env.example` is committed; `.env` is git-ignored. In production, config comes
from AWS Secrets Manager (chapter 11) — never from a file in the image.

## Where a bug gets fixed

The most useful maintenance habit in this codebase. When something is wrong, ask
**which layer owns this rule?**

| Symptom | Fix belongs in |
|---|---|
| Year 3000 is accepted | `server/src/domain/Year.js` |
| Empty country list where events exist | `server/src/services/TimelineService.js` |
| Wrong HTTP status | `server/src/api/errorHandler.js` or the controller |
| Wrong sort order | the repositories — **and the contract test** |
| Events shown for the wrong year | `client/js/App.js` state |
| Focus lost while navigating the timeline | `client/js/components/TimelineView.js` |
| Text renders as markup | `client/js/dom.js` — and treat it as a security incident |
| An uncited event got in | `server/src/domain/HistoricalEvent.js` — and ask how it got past |

Fixing a symptom in the wrong layer is how a codebase rots: the rule ends up
duplicated in three places, and the next person changes two of them.

## Code review

What to look for, roughly in order:

1. **Does it do what it says?** Read the tests first — they state the intent.
2. **Is the rule in the right layer?** (The table above.)
3. **What happens on the unhappy path?** Empty, missing, malformed, slow.
4. **Is it testable?** If it is hard to test, it is usually too coupled.
5. **Will the next person understand it?** Comments explain *why*; the code says
   *what*.

Two review questions specific to this stack, both with yes/no answers:

```bash
grep -rn '\${' server/src/repositories/postgres/*.js   # interpolation in SQL?
grep -rn 'innerHTML' client/js/                        # unescaped markup?
```

Today the first returns two `${EVENT_COLUMNS}` hits (a hard-coded constant, not
user input) and the second returns only the warning comments in `dom.js`. Any
*new* hit in either is a review conversation, not a nit.

Not worth a comment: formatting, or style preferences with no behavioural
difference. Reviews that bikeshed on naming while missing an N+1 query are worse
than no review.

## Documentation that stays true

Documentation lies as soon as the code changes. The defence is to **make as much
of it executable as possible**:

- Tests document behaviour and fail when it changes.
- Comments explain *why* — the part that does not go stale, because reasons
  outlive implementations.

```js
// ❌ says what the next line says, and will outlive its truth
// Set the year
this.year = year;

// ✅ explains a decision the code cannot
// COUNT(*) returns BIGINT, which exceeds JavaScript's safe integer range, so
// `pg` hands it back as a STRING rather than silently losing precision.
```

## Performance, when it is time

**Measure first.** The three things most likely to matter here:

- **N+1 queries** — a query inside a loop. `countByYearRange` exists to avoid one.
- **Missing indexes** — `EXPLAIN ANALYZE` on your slowest query.
- **Unnecessary DOM work** — the reason `TimelineView` updates two elements
  instead of 126.

Then caching: `/api/timeline` changes rarely and would happily sit behind
`Cache-Control: public, max-age=300` and a CDN.

Do not optimise before measuring. The bottleneck is rarely where you think.

## Do it yourself

1. **Read the CI file end to end.** Explain what each step protects against.
2. **Make CI fail.** Push a branch with a deliberately broken test.
3. **Add ESLint** and fix what it finds.
4. **Add request logging** middleware: method, path, status, duration, in JSON.
5. **Audit the dependencies.** `npm outdated`, then update one minor version and
   run the suite.

## Check yourself

- Why does CI run the Postgres contract tests when local `npm test` skips them?
- Why validate configuration at startup rather than on first use?
- What replaces a compiler in a plain-JavaScript project?

<details>
<summary>Answers</summary>

- **CI runs them**: the local loop optimises for speed (no database required); CI
  optimises for confidence and has a Postgres service container. The skip is a
  local convenience, not a permanent hole.
- **Startup validation**: a config error at boot is one obvious failure before any
  traffic arrives. The same error on first use is an outage at 3am, reported as
  "the site is broken sometimes".
- **Replacing a compiler**: validation at the boundaries, contract tests for every
  abstraction with more than one implementation, a linter, and — when the project
  outgrows what tests comfortably cover — TypeScript.

</details>

→ Next: [11 — Deploying to AWS](11-aws-deployment.md)
