# 02 — TypeScript and OOP

> **Concept** → The code → Do it yourself → Check yourself

## Why TypeScript at all

You asked for OOP as it is used in industry. That decided the language: classes,
interfaces and dependency inversion only really *teach* with types, because
without them "this class implements that contract" is a comment rather than a
fact the compiler checks.

TypeScript is JavaScript plus a type layer that is erased before the code runs.
Nothing in `dist/` has types in it. The types exist to catch mistakes at the one
moment they are cheapest to fix: while you are typing.

```ts
// JavaScript: this is fine until 3am, when a user hits it.
function eventsFor(year, country) { ... }
eventsFor('FRA', 1969);   // arguments swapped. Runs. Returns nothing. No error.

// TypeScript: this never reaches a user.
function eventsFor(year: number, country: string) { ... }
eventsFor('FRA', 1969);
//        ~~~~~ Argument of type 'string' is not assignable to parameter of type 'number'.
```

## The four pillars, in this codebase

Every OOP course lists these four. Here is where each one actually earns its
keep in this project — not as a definition, but as a file you can open.

### 1. Encapsulation — hide the data, expose the rules

`src/domain/Year.ts`. A `number` can be `-4`, `NaN`, or `99999`. A `Year` cannot:

```ts
export class Year {
  readonly value: number;

  constructor(value: number) {
    if (!Number.isInteger(value)) throw new ValidationError(...);
    if (value < Year.EARLIEST)    throw new ValidationError(...);
    if (value > Year.latestAllowed()) throw new ValidationError(...);
    this.value = value;
    Object.freeze(this);
  }
}
```

The payoff: **once a `Year` exists, every function downstream can trust it
without re-checking.** The validation happens once, at the boundary. This idea
has a name — *parse, don't validate* — and it is the single most useful habit in
this chapter.

`Config` in `src/config/Config.ts` takes it further with a **private
constructor**: the only way to get a `Config` is `Config.fromEnv()`, so an
unvalidated one is not merely discouraged, it is unconstructable.

### 2. Inheritance — share behaviour, but sparingly

`src/api/errorHandler.ts` handles every domain error with one branch:

```ts
if (error instanceof DomainError) {
  res.status(error.httpStatus).json({ error: { code: error.code, ... } });
}
```

`ValidationError` and `NotFoundError` both extend `DomainError`. Adding a
`ConflictError` with `httpStatus = 409` requires **zero changes** to the handler.

On the client, `TimelineApi extends ApiClient` so that fetch, JSON parsing,
error mapping and abort handling are written once and every API method is a
one-liner.

> **The warning that belongs with every inheritance lesson:** inheritance is the
> tightest coupling in object-oriented programming. A subclass depends on its
> parent's internals, and the dependency is invisible at the call site. Both
> hierarchies here are **one level deep and used for genuine "is-a"
> relationships**. When you find yourself three levels deep, you almost always
> wanted composition instead — which is why `InMemoryDataset` is *given* to the
> repositories rather than inherited by them.

### 3. Polymorphism — one call site, many implementations

`TimelineService` calls `this.countries.findAll()`. At runtime that might be
Postgres or two arrays. The service cannot tell and does not care.

This is what makes the service testable, and it is the mechanism behind the
contract test in chapter 5.

### 4. Abstraction — depend on the shape, not the thing

```ts
export class TimelineService {
  constructor(
    private readonly countries: CountryRepository,   // abstract
    private readonly events: EventRepository,        // abstract
  ) {}
}
```

Never `PostgresCountryRepository`. The concrete classes are named in exactly one
file: `src/main.ts`.

## `interface` or `abstract class`?

TypeScript gives you both. The distinction matters and is worth getting right:

| | `interface` | `abstract class` |
|---|---|---|
| Exists at runtime | ❌ erased | ✅ real JS class |
| `instanceof` works | ❌ | ✅ |
| Can hold shared code | ❌ | ✅ |
| Multiple inheritance | ✅ | ❌ single only |

This project uses:

- **`abstract class`** for repositories (`CountryRepository`). They are extended,
  they benefit from `instanceof`, and they may gain shared helpers later.
- **`interface`** for plain data shapes (`CountryProps`, `TimelineTick`,
  `QueryRunner`). No behaviour, no runtime cost.

Rule of thumb: **interface for data, abstract class for something you extend.**

## SOLID, one file each

Not as an acronym to recite — as five specific files in this repository.

**S — Single responsibility.** `TimelineService` knows what the product does.
`PostgresEventRepository` knows SQL. `TimelineController` knows HTTP. Ask "what
would make me edit this file?" If there are two unrelated answers, split it.

**O — Open/closed.** Open to extension, closed to modification. Add
`ConflictError extends DomainError` and the error handler supports it without
being edited.

**L — Liskov substitution.** Any `CountryRepository` must work anywhere a
`CountryRepository` is expected. This is easy to *say* and easy to *break* — one
implementation sorts, the other does not, and something breaks in production
only. `tests/unit/repositoryContract.test.ts` runs the same assertions against
both implementations to keep them honest. **A shared contract test is what turns
Liskov from a slogan into something enforced.**

**I — Interface segregation.** `CountryRepository` has four methods, all about
countries. A fat `DataRepository` with forty methods would force
`InMemoryCountryRepository` to implement things it never uses.

**D — Dependency inversion.** Both `TimelineService` and
`PostgresCountryRepository` depend on abstractions (`CountryRepository`,
`QueryRunner`), not on each other's concrete types.

## Where this project does *not* use classes

Honesty matters more than consistency: **React components here are functions,
not classes.** React moved to function components and hooks years ago, class
components are effectively legacy, and writing them would teach you something
you would have to unlearn.

So the rule this project actually follows is:

> **Classes for things with identity and behaviour** — a `Year`, an
> `ApiClient`, a repository. **Functions for transformations** — a component
> that maps props to markup, `formatMonthDay`, a reducer.

That is how modern TypeScript codebases are written. "OOP as industry standard"
means using objects where objects fit, not making everything a class.

## TypeScript settings worth knowing

From `server/tsconfig.json` — each of these has caught a real bug in this
project:

```jsonc
"strict": true,                  // the whole strict family. Never turn this off.
"noUncheckedIndexedAccess": true // arr[0] is T | undefined — because it is
"noImplicitOverride": true,      // must write `override` — catches renamed methods
```

`noUncheckedIndexedAccess` is the one people disable first and regret. It is why
`src/scripts/seed.ts` reads `rows[0]?.count ?? '0'` instead of `rows[0].count`.
That query returns a row *today*; the type says "prove it".

## Do it yourself

1. **Make the compiler catch a real bug.** In `TimelineService`, change
   `getEventsForYearAndCountry(year: number, countryCode: string)` to
   `(countryCode: string, year: number)` — swap them — but do not touch the
   controller. Run `npx tsc --noEmit`. Read the error. Now imagine that same
   change in untyped JavaScript.

2. **Add an error type.** Create `ConflictError extends DomainError` with
   `code = 'CONFLICT'` and `httpStatus = 409`. Throw it from a service method.
   Confirm the API returns 409 **without editing `errorHandler.ts`**. That is
   open/closed, demonstrated rather than defined.

3. **Break Liskov deliberately.** Delete the `.sort()` from
   `InMemoryCountryRepository.findAll()`. Run `npx vitest run
   repositoryContract`. Watch the contract test catch it. This is the exercise
   that makes chapter 5 land.

## Check yourself

- Why is `Year` immutable, and what would break if it were not?
- Why does `Country` compare by `id` while `Year` compares by `value`?
- Why is `Config`'s constructor private?
- When would you choose composition over inheritance?

<details>
<summary>Answers</summary>

- **Immutability**: a `Year` can be shared, cached and used as a map key
  safely. If any holder could mutate it, a `Year` passed into a function could
  come back different — and `Object.freeze` means that fails loudly instead of
  silently.
- **`Country` vs `Year`**: `Country` is an *entity* — it has identity that
  survives its attributes changing (fix a typo in the name, still the same
  country). `Year` is a *value object* — 1969 is 1969; there is no "which 1969".
  Chapter 3 is entirely about this distinction.
- **Private constructor**: it makes `Config.fromEnv()` the only way in, so an
  unvalidated `Config` cannot exist anywhere in the program.
- **Composition over inheritance**: when the relationship is "has-a" or "uses-a"
  rather than "is-a", and whenever you would otherwise inherit just to reuse a
  method. `InMemoryCountryRepository` *has* a dataset; it is not *a* dataset.

</details>

→ Next: [03 — The domain model](03-domain-model.md)
