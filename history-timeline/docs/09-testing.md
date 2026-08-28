# 09 — Testing

> **Concept** → The code → Do it yourself → Check yourself

## What this project actually has

```
server: 89 tests  (+ 9 opt-in Postgres contract tests)
client: 74 tests
───────────────────────────────────────────────────────
        163 tests in about 1.5 seconds
```

Speed is not vanity. **A suite you can run on every save is a different tool from
one you run before lunch.** Under two seconds and you run it constantly; over
thirty and you stop.

## No test framework

```bash
node --test server/tests/*.test.js client/tests/*.test.js
```

That is the whole setup. Node has had a built-in test runner since v18:

```js
import { describe, it, before, beforeEach } from 'node:test';
import assert from 'node:assert/strict';
```

No Jest, no Vitest, no config file. One fewer tool to learn while you are
learning everything else — and it makes the point that **testing is not magic: a
test is a function that throws when something is wrong.**

`assert/strict` (rather than plain `assert`) uses `===` semantics for
`deepEqual`, so `'1'` never passes as `1`.

The one thing we do install is **jsdom**, because Node has no `document`. That is
the entire trick behind every front-end test runner; `client/tests/domEnvironment.js`
does it in fifteen visible lines rather than in a config file.

## Testing matters more without a compiler

In a typed language, the compiler catches whole categories of mistakes for free:
a renamed method, a swapped argument, a missing implementation.

In plain JavaScript, **nothing does that but your tests.** So the balance shifts:
tests here are not only checking behaviour, they are also doing the job a type
checker would. The clearest example is the repository contract test — a compiler
would notice `findByRegion` was missing from one implementation; here, the test
is what notices.

**The looser the language, the more the tests have to carry.** That is a real
trade-off of this stack, and worth naming honestly.

## The pyramid, as built here

```
        ╱ App.test.js  ╲            few, slower, high confidence
       ╱  (whole app)   ╲
      ╱──────────────────╲
     ╱  api.test.js       ╲         integration: full HTTP stack
    ╱   views.test.js      ╲
   ╱────────────────────────╲
  ╱ domain, Store, dom,       ╲     many, fast, precise
 ╱  TimelineService, ApiClient  ╲
╱────────────────────────────────╲
```

| Level | Count | What it proves | Speed |
|---|---|---|---|
| Unit | ~120 | One rule in isolation | µs |
| Integration | ~36 | Layers fit together | ms |
| Whole-app | 7 | The features work end to end | ~100 ms |

**Why more unit tests?** When a unit test fails you know exactly which rule
broke. When a whole-app test fails you know *something* broke. Both are useful;
the precise ones should be plentiful.

## What to test — and what never to

**Test:** business rules, boundaries, error paths, contracts between layers, and
every bug you fix (a test is how a bug stays fixed).

**Do not test:** implementation details. The heuristic — **would this test still
pass if I rewrote the internals but kept the behaviour?** If no, you are testing
the wrong thing.

```js
// ❌ implementation — breaks on any refactor
assert.equal(view._buttons.length, 4);

// ✅ behaviour — survives any rewrite that keeps the feature
assert.equal(root.querySelectorAll('[role="radio"]').length, 4);
```

There is one deliberate exception in this project, and it is worth understanding
because it shows the rule is a heuristic rather than a law:

```js
it('reuses the same button elements when only the selection changes', () => {
  const before = root.querySelector('[data-year="1968"]');
  view.update(ready(1968));
  assert.equal(before, root.querySelector('[data-year="1968"]'));
});
```

That asserts on DOM node identity — an implementation detail by any normal
reading. But **the behaviour it protects is real** (keyboard focus and scroll
position surviving an update), and it is otherwise only observable in a real
browser. When an implementation detail is the only observable proxy for a
behaviour you care about, testing it is right — but say so in a comment, as that
file does.

## Boundary testing

Bugs live at edges. `domain.test.js` tests every one:

```js
it('accepts the earliest supported year, 100 AD', ...);   // the edge
it('rejects years before 100 AD', ...);                   // edge − 1
it('accepts the current year', ...);                      // the other edge
it('rejects future years', ...);                          // edge + 1
```

For any range: below, at, just inside, just below the top, at the top, above.
That is where off-by-one lives.

Without `it.each` (Node's runner has no such helper), a loop does the job — and
the message argument keeps failures readable:

```js
for (const bad of [1969.5, Number.NaN, '1969', null, undefined, {}]) {
  assert.throws(() => new Year(bad), ValidationError, `should reject ${JSON.stringify(bad)}`);
}
```

**Always pass that third argument.** Without it, a failure says "expected to
throw" and you have to work out which of six inputs it was.

## Test doubles: the honest hierarchy

Ranked by how much confidence they give:

1. **The real thing** — best, when it is fast. The domain tests use real `Year`
   and `HistoricalEvent` objects.
2. **A real alternative implementation** — `InMemoryCountryRepository` is real
   code satisfying the real contract, verified by the contract test.
3. **A stub at the system boundary** — `stubFetch` replaces the network only.
   `ApiClient`'s real parsing and error handling still run.
4. **A mock of your own class** — last resort. Replacing `TimelineApi` would mean
   the tests never execute the code you ship.

This project uses 1, 2 and 3. **The further down you go, the more you are testing
your assumptions instead of your code.**

Concretely, on the client we stub `globalThis.fetch`, so a bug in `ApiClient`'s
502-handling is caught:

```js
it('survives an error response that is not JSON', async () => {
  active = stubFetch({ '/api/timeline': '<html>502 Bad Gateway</html>' },
                     { status: { '/api/timeline': 502 } });
  const error = await new TimelineApi().getTimeline().catch((e) => e);
  assert.equal(error.code, 'UNKNOWN_ERROR');
  assert.equal(error.isClientError, false);
});
```

A proxy returning HTML is a real thing that happens in production.

Note the stub returns a `restore()` function and every test calls it. **Leaving a
global replaced leaks into the next test**, and those failures are miserable to
debug.

## HTTP tests without a library

No supertest. Start the real app on port 0 — "any free port" — and use Node's
built-in `fetch`:

```js
before(async () => {
  const app = createApp({ service, db: null, serveClient: false });
  await new Promise((resolve) => { server = app.listen(0, resolve); });
  baseUrl = `http://127.0.0.1:${server.address().port}`;
});
```

Fewer dependencies, and it makes clear that an HTTP test is just an HTTP request.

`serveClient: false` matters: with the static handler on, an unknown path would
fall through to it instead of returning the JSON 404 the test checks.

## Contract tests

Two implementations of one abstraction: write the assertions **once**, run them
against **both**. Fast feedback by default, real confidence on demand:

```bash
npm run db:up && npm run db:migrate && npm run db:seed
TEST_DATABASE_URL=postgres://historian:historian@localhost:5432/history_timeline npm test
```

Node's runner takes a `skip` option, so the Postgres block announces itself as
skipped rather than silently vanishing:

```js
describe('PostgreSQL repositories', { skip: !testDatabaseUrl && 'set TEST_DATABASE_URL to run' }, ...)
```

## Testing the *content*

For an educational product the data **is** the product, so `seedData.test.js`
protects the editorial rules:

```js
it('every event is valid, cited, and in range', ...);
it('every event belongs to a known country', ...);
it('no country has a duplicate event in the same year', ...);
it('every country has at least three events, so no country is a dead end', ...);
```

Add an event with a broken URL and the build fails — before a learner sees it.

**Whatever your product's equivalent is, test it.** The rules that define
correctness are not always in the code.

## Testing the DOM as a user

```js
[...root.querySelectorAll('button')].find((b) => b.textContent.includes('Japan')).click();
assert.deepEqual(chosen, ['JPN']);
```

Queries by **role** and **visible text** — the way a user (or a screen reader)
finds things.

A side effect worth noticing: **querying by role forces accessible markup.** If
`querySelectorAll('[role="radio"]')` finds nothing, that is not the test being
awkward — the test just told you a keyboard user cannot reach it. The test suite
and the accessibility audit are the same activity.

### The security test

```js
it('renders a malicious title as text, never as markup', () => {
  const attack = '<img src=x onerror="globalThis.__pwnedView = true">';
  view.update(state({ events: Async.success({ ...EVENTS, events: [anEvent({ title: attack })] }) }));
  assert.equal(root.querySelector('img'), null);
  assert.equal(globalThis.__pwnedView, undefined);
});
```

Keep this kind of test forever. The unsafe version **looks fine** until the day
someone exploits it, so a test is the only thing standing between a refactor and
a vulnerability.

## Async tests, and the `flush` helper

The app loads data in `async` methods, so after a click the DOM is not updated
until those promises settle:

```js
button('France').click();
await flush();                       // let pending promises run
assert.match(text('events'), /Concorde first flight/);
```

`flush` yields to the event loop a few times. It is enough for promise chains
that are not waiting on real I/O — and since the network is stubbed, there is
none.

**Never use a fixed `setTimeout(500)` to "wait for the UI".** It is slow when it
works and flaky when it does not.

## Coverage, used correctly

```bash
npm run test:coverage
```

**Coverage tells you what is definitely untested. It does not tell you what is
well tested.** 100% coverage with no assertions is 100% meaningless.

Use it to find blind spots — "the whole error path is uncovered" is a real
finding. Do not chase a number, and never make it a gate that people game.

## Naming tests

A test name should read as a sentence about the product:

```js
✅ it('refuses to exist without a citable source')
✅ it('returns an empty list — not an error — for a year with no events')
❌ it('test event validation')
❌ it('works')
```

When a failure appears in CI at 6pm, the name is all you have.

## Do it yourself

1. **Break things and read the failures.** Change `Year.EARLIEST` to 1500. Which
   tests fail, and are the messages good enough to diagnose from?

2. **Test-drive a feature.** Write a failing test for "events can be filtered by
   category" *first*. Watch it fail. Implement it. Watch it pass. Doing that loop
   once on something real is worth more than reading about TDD.

3. **Fix a bug properly.** Find any bug, write the test that reproduces it,
   *then* fix it. Confirm the test fails before and passes after.

4. **Run the contract tests against Postgres** (commands above) and confirm all
   18 pass.

5. **Delete a repository method** from `InMemoryCountryRepository` and run the
   tests. Notice the base class's "must implement" error naming the exact class
   and method — that message is doing a compiler's job.

## Check yourself

- Why stub `fetch` rather than replace `TimelineApi`?
- Why are the Postgres contract tests opt-in rather than always on?
- What makes a test brittle?

<details>
<summary>Answers</summary>

- **Stub at the boundary**: `ApiClient`'s real error handling, JSON parsing and
  abort logic then run in the tests. Replacing `TimelineApi` would mean the tests
  never execute the code you actually ship.
- **Opt-in**: they need a running, migrated, seeded database. Requiring that for
  `npm test` would make the everyday loop slow and fragile. CI sets
  `TEST_DATABASE_URL`, so nothing is skipped where it counts.
- **Brittle**: it asserts on *how* rather than *what* — a private field, a call
  order, an exact string free to change. It fails on refactors that broke
  nothing and teaches the team to distrust the suite.

</details>

→ Next: [10 — Maintaining it](10-maintenance.md)
