# 08 — Testing

> **Concept** → The code → Do it yourself → Check yourself

## What this project actually has

```
server: 97 tests + 9 opt-in Postgres contract tests   ~1.0 s
client: 35 tests                                       ~1.9 s
────────────────────────────────────────────────────────────
        132 tests in about three seconds
```

Speed is not vanity. **A suite you can run on every save is a different tool
from one you run before lunch.** Under two seconds and you run it constantly;
over thirty and you stop.

## The pyramid, as built here

```
        ╱ App.test.tsx ╲            few, slow, high confidence
       ╱  (whole app)   ╲
      ╱──────────────────╲
     ╱  api.test.tsx      ╲         integration: full HTTP stack
    ╱   components.test    ╲
   ╱────────────────────────╲
  ╱  Year, entities,          ╲     many, fast, precise
 ╱   TimelineService, ApiClient ╲
╱────────────────────────────────╲
```

| Level | Count | What it proves | Speed |
|---|---|---|---|
| Unit | ~90 | One rule in isolation | µs |
| Integration | ~30 | Layers fit together | ms |
| Whole-app | 6 | The features work end to end | ~100 ms |

**Why more unit tests?** When a unit test fails you know exactly which rule
broke. When a whole-app test fails you know *something* broke. Both are useful;
the precise ones should be plentiful.

## What to test — and what never to

**Test:** business rules, boundaries, error paths, contracts between layers, and
every bug you fix (a test is how a bug stays fixed).

**Do not test:** implementation details. Nothing in this project asserts on a CSS
class name, a private method, or a component's internal state. Those tests break
on every refactor while catching none of the bugs that matter.

The heuristic: **would this test still pass if I rewrote the internals but kept
the behaviour?** If no, you are testing the wrong thing.

Compare:

```ts
// ❌ implementation — breaks the moment you rename a class
expect(wrapper.find('.timeline__year--selected')).toHaveLength(1);

// ✅ behaviour — survives any refactor that keeps the feature
expect(screen.getByRole('radio', { checked: true })).toHaveAccessibleName(/1969/);
```

## Boundary testing

Bugs live at edges. `Year.test.ts` tests every one:

```ts
it('accepts the earliest supported year, 100 AD', ...);   // the edge
it('rejects years before 100 AD', ...);                   // edge − 1
it('accepts the current year', ...);                      // the other edge
it('rejects future years', ...);                          // edge + 1
```

For any range, test: below, at, just inside, just below the top, at the top,
above. That is where off-by-one lives.

`it.each` keeps this from becoming tedious:

```ts
it.each(['', 'nineteen-sixty-nine', '19a9', '1969.0'])('rejects %o', (raw) => {
  expect(() => Year.fromString(raw)).toThrow(ValidationError);
});
```

Four tests, four separate failure messages, four lines.

## Test doubles: the honest hierarchy

Ranked by how much confidence they give:

1. **The real thing** — best, when it is fast. The domain tests use real `Year`
   and `HistoricalEvent` objects.
2. **A real alternative implementation** — `InMemoryCountryRepository` is real
   code that satisfies the real contract, verified by the contract test.
3. **A stub at the system boundary** — `stubFetch` replaces the network only.
   `ApiClient`'s real parsing and error handling still run.
4. **A mock of your own class** — last resort. `vi.mock('./TimelineApi')` would
   mean the tests never execute the code you ship.

This project uses 1, 2 and 3. **The further down you go, the more you are testing
your assumptions instead of your code.**

Concretely, on the client we stub `globalThis.fetch` rather than mocking
`TimelineApi`, so a bug in `ApiClient`'s 502-handling is caught:

```ts
it('survives an error response that is not JSON', async () => {
  vi.spyOn(globalThis, 'fetch').mockResolvedValue(
    new Response('<html>502 Bad Gateway</html>', { status: 502 }),
  );
  const error = await new TimelineApi('').getTimeline().catch((e) => e) as ApiError;
  expect(error.code).toBe('UNKNOWN_ERROR');
  expect(error.isClientError).toBe(false);
});
```

A proxy returning HTML is a real thing that happens in production.

## Contract tests (chapter 5, from the testing angle)

Two implementations of one abstraction: write the assertions **once**, run them
against **both**.

```ts
contractFor('InMemory repositories', ...);                 // always
describe.skipIf(!testDatabaseUrl)('PostgreSQL repositories', () => {
  contractFor('Postgres repositories', ...);               // opt-in
});
```

Fast feedback by default; real confidence on demand:

```bash
npm run db:up && npm run db:migrate && npm run db:seed
TEST_DATABASE_URL=postgres://historian:historian@localhost:5432/history_timeline npm test
```

## Testing the *content*

For an educational product the data **is** the product, so `seedData.test.ts`
protects the editorial rules:

```ts
it('every event is valid, cited, and in range', ...);
it('every event belongs to a known country', ...);
it('no country has a duplicate event in the same year', ...);
it('every country has at least three events, so no country is a dead end', ...);
```

Add an event with a broken URL and the build fails — before a learner sees it.

**Whatever your product's equivalent is, test it.** The rules that define
correctness are not always in the code.

## Testing components as a user

```ts
await userEvent.click(screen.getByRole('button', { name: /France/ }));
expect(await screen.findByText('Concorde first flight')).toBeInTheDocument();
```

Queries by **role** and **visible text** — the way a user (or a screen reader)
finds things.

A side effect worth noticing: **querying by role forces accessible markup.** If
`getByRole('button')` cannot find your clickable `div`, that is not the test
being awkward — the test just told you a keyboard user cannot reach it. The
test suite and the accessibility audit are the same activity.

`findBy*` waits for async updates; `getBy*` throws immediately; `queryBy*`
returns null and is how you assert something is *absent*:

```ts
await waitFor(() => {
  expect(screen.queryByText('Concorde first flight')).not.toBeInTheDocument();
});
```

## The test that found a real bug

Worth repeating because it is the best argument for writing tests at all:

```ts
it('clears the chosen country when the year changes', async () => {
  await userEvent.click(await screen.findByRole('button', { name: /France/ }));
  expect(await screen.findByText('Concorde first flight')).toBeInTheDocument();

  await userEvent.click(screen.getByRole('radio', { name: /1967/ }));

  await waitFor(() => {
    expect(screen.queryByText('Concorde first flight')).not.toBeInTheDocument();
  });
});
```

This failed on first run. The cause: a disabled `useAsyncResource` kept its last
`success` state, so France's 1969 events stayed on screen under a 1967 heading —
**plausible-looking, completely wrong data, and invisible by eye.**

The fix was an `idle` state in the hook (chapter 7). The test was written from
the *requirement* ("changing the year clears the country"), not from the
implementation — which is exactly why it could catch a bug the implementation
did not anticipate.

## Coverage, used correctly

```bash
npm run test:coverage --workspace=server
```

**Coverage tells you what is definitely untested. It does not tell you what is
well tested.** 100% coverage with no assertions is 100% meaningless.

Use it to find blind spots — "the whole error path is uncovered" is a real
finding. Do not chase a number, and never make it a gate that people game.

## Naming tests

A test name should read as a sentence about the product:

```ts
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
   category" *first*. Watch it fail. Implement it. Watch it pass. That loop is
   TDD, and doing it once on something real is worth more than reading about it.

3. **Fix a bug properly.** Find any bug, write the test that reproduces it,
   *then* fix it. Confirm the test fails before and passes after. This is the
   habit that makes a suite valuable over years.

4. **Run the contract tests against Postgres** (commands above) and confirm all
   18 pass.

## Check yourself

- Why stub `fetch` rather than mock `TimelineApi`?
- Why are the Postgres contract tests opt-in rather than always on?
- What makes a test brittle?

<details>
<summary>Answers</summary>

- **Stub at the boundary**: `ApiClient`'s real error handling, JSON parsing and
  abort logic then run in the tests. Mocking `TimelineApi` would mean the tests
  never execute the code you actually ship.
- **Opt-in**: they need a running, migrated, seeded database. Requiring that for
  `npm test` would make the everyday loop slow and fragile. CI runs them with a
  Postgres service container (chapter 9), so nothing is skipped where it counts.
- **Brittle**: it asserts on *how* rather than *what* — a class name, a private
  method, a call order, an exact string that is free to change. It fails on
  refactors that broke nothing and, worse, teaches the team to distrust the suite.

</details>

→ Next: [09 — Maintaining it](09-maintenance.md)
