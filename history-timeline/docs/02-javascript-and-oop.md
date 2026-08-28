# 02 — JavaScript and OOP

> **Concept** → The code → Do it yourself → Check yourself

This chapter is the JavaScript you need for the rest of the course, taught
through the code in this project rather than in the abstract.

## Modules: `import` and `export`

Every file in this project is an **ES module**. Two things make that work:

```json
// package.json
{ "type": "module" }
```

```html
<script type="module" src="/js/app.js"></script>
```

```js
export class Year { ... }              // named export
import { Year } from './Year.js';      // named import
```

Three rules that trip people up:

1. **The `.js` extension is required.** `from './Year'` fails in the browser.
   Node and browsers resolve real file paths; the extensionless style you see
   elsewhere is a bundler feature, and we have no bundler.
2. **Modules are always in strict mode.** Assigning to an undeclared variable
   throws instead of silently creating a global.
3. **Top-level names are module-scoped**, not global. Two files can both have a
   `const config` without colliding.

You will also meet the older **CommonJS** style (`require`, `module.exports`) in
tutorials and older packages. It is the previous system. This project uses ES
modules throughout because they work in the browser too, unchanged.

## Classes

```js
export class Year {
  #value;                        // private field
  static EARLIEST = 100;         // static (class-level) field

  constructor(value) { ... }     // runs on `new Year(1969)`

  get value() { return this.#value; }        // getter: year.value, not year.value()
  static fromString(raw) { ... }             // called on the class: Year.fromString('1969')
  equals(other) { ... }                      // instance method: year.equals(other)
}
```

`extends` and `super`:

```js
export class ValidationError extends DomainError {
  constructor(message, field) {
    super(message);              // MUST call super() before using `this`
    this.field = field;
  }
}
```

## `#private` is real privacy

This is the modern feature that makes JavaScript encapsulation genuine rather
than a naming convention:

```js
class Year {
  #value;
  get value() { return this.#value; }
}

const y = new Year(1969);
y.value        // 1969  — through the getter
y.#value       // SyntaxError: Private field '#value' must be declared in an enclosing class
Object.keys(y) // []    — invisible to iteration and JSON.stringify
```

Compare the older convention, `this._value`, which is just a name — anyone can
write `y._value = 3050` and nothing stops them.

**Both privacy techniques appear in this project on purpose.** `Year` uses
`#value` plus a getter. `Country` and `HistoricalEvent` use public fields plus
`Object.freeze(this)`. Freezing is weaker (fields are readable, and it is
shallow) but far less boilerplate for a four-field object, and it delivers what
we need: nobody can change it after construction. Choosing the lighter tool when
it is sufficient is a real engineering skill.

## `this`, and the bug you will definitely hit

`this` in JavaScript is decided by **how a function is called**, not where it was
written. Detach a method from its object and `this` is gone:

```js
class TimelineController {
  getTimeline(req, res) { this.#service...; }   // a normal method
}

router.get('/timeline', controller.getTimeline);
//                      ^ passed as a plain function — `this` is undefined inside
```

Three fixes; this project uses the third:

```js
router.get('/timeline', (req, res) => controller.getTimeline(req, res)); // wrap
constructor() { this.getTimeline = this.getTimeline.bind(this); }        // bind
#getTimeline = async (req, res) => { ... };                              // arrow field ✅
```

An **arrow function has no `this` of its own** — it uses the `this` of the scope
where it was written. As a class field, that scope is the constructor, so the
reference is safe to pass anywhere. Same technique in `TimelineView`, where
handlers are passed to `addEventListener`.

This is the single most common "why is `this` undefined" bug in JavaScript.
Recognising it will save you hours.

## Runtime validation is your type system

Without a compiler, nothing catches this:

```js
getEventsForYearAndCountry('FRA', 1969)   // arguments swapped
```

So the constructor does:

```js
constructor(value) {
  if (typeof value !== 'number' || !Number.isInteger(value)) {
    throw new ValidationError(`Year must be a whole number, received ${JSON.stringify(value)}.`);
  }
  ...
}
```

**Every rule you do not write is a rule that does not exist.** That is the
central discipline of writing plain JavaScript well, and it is why the domain
layer in chapter 3 is as strict as it is.

You will notice `typeof value !== 'number'` as well as `Number.isInteger`.
`Number.isInteger('1969')` is already `false`, but the explicit `typeof` check
makes the intent obvious and the error message accurate.

> **When you later add TypeScript**, this is the natural next step and the code
> is ready for it — the constructors stay, because a type checker validates what
> your *code* does, not what a *user* sends. Validation at the boundary is
> needed in both worlds.

## The four pillars, in this codebase

Not as definitions — as files you can open.

### Encapsulation — hide the data, expose the rules

`src/domain/Year.js`. A `number` can be `-4`, `NaN` or `"1969"`. A `Year` cannot.
Once one exists, every function downstream can trust it. That idea is called
*parse, don't validate*, and it is the most useful habit in this chapter.

`Config` goes further with a **private constructor** (`src/config/Config.js`):
`Config.fromEnv()` is the only way in, so an unvalidated config cannot exist.

### Inheritance — share behaviour, sparingly

`src/api/errorHandler.js` handles every domain error in one branch:

```js
if (error instanceof DomainError) {
  res.status(error.httpStatus).json({ error: { code: error.code, ... } });
}
```

Adding `ConflictError extends DomainError` with `httpStatus = 409` requires
**zero changes** here.

> **The warning that belongs with every inheritance lesson:** inheritance is the
> tightest coupling in OOP, and the dependency is invisible at the call site.
> Both hierarchies here are one level deep and genuine "is-a" relationships.
> Three levels deep, you almost always wanted composition — which is why
> `InMemoryDataset` is *given* to the repositories rather than inherited.

### Polymorphism — one call site, many implementations

`TimelineService` calls `this.#countries.findAll()`. At runtime that might be
Postgres or two arrays. The service cannot tell and does not care.

### Abstraction — depend on the shape, not the thing

```js
constructor(countryRepository, eventRepository) { ... }  // abstract types
```

Never `PostgresCountryRepository`. Concrete classes are named in exactly one
file: `src/main.js`.

## Expressing an interface without a compiler

JavaScript has no `interface` and no `abstract`. Here is how this project says
both anyway:

```js
export class CountryRepository {
  constructor() {
    if (new.target === CountryRepository) {
      throw new TypeError('CountryRepository is abstract; extend it.');
    }
  }
  async findAll() {
    throw new Error(`${this.constructor.name} must implement findAll()`);
  }
}
```

`new.target` is the constructor that was actually called — `CountryRepository`
for a direct `new`, `InMemoryCountryRepository` for a subclass. The same trick
gives every `DomainError` its correct `name`.

**The honest limitation:** a typed language catches a missing method at compile
time; this catches it the first time the method is called. That is exactly why
`tests/repositoryContract.test.js` exists — it calls every method on every
implementation, so "the first time it is called" happens in CI, not in
production. **The looser the language, the more the tests have to carry.**

## SOLID, one file each

**S — Single responsibility.** `TimelineService` knows the product.
`PostgresEventRepository` knows SQL. `TimelineController` knows HTTP. Ask "what
would make me edit this file?" Two unrelated answers means split it.

**O — Open/closed.** Add a `DomainError` subclass; the error handler supports it
unedited.

**L — Liskov substitution.** Any `CountryRepository` must work anywhere one is
expected. Easy to say, easy to break — one implementation sorts, the other does
not. `tests/repositoryContract.test.js` runs the same assertions against both.
**A shared contract test is what turns Liskov from a slogan into something
enforced.**

**I — Interface segregation.** `CountryRepository` has four methods, all about
countries. `HealthController` depends on "something with `isHealthy()`", not on
a whole database class.

**D — Dependency inversion.** Services depend on abstract repositories;
repositories depend on a `query`-shaped object, not on `pg` directly.

## Where this project does *not* use classes

Honesty over consistency. Classes are used for **things with identity and
behaviour** — a `Year`, an `ApiClient`, a view, a repository. Plain functions are
used for **transformations**:

```js
export function formatMonthDay(monthDay) { ... }   // pure: input → output
export function groupByRegion(countries) { ... }   // pure
export function el(tag, props, children) { ... }   // a factory, not a class
```

Pure functions are the easiest thing in software to test, which is why those are
exported separately rather than buried inside a class. "OOP as industry
standard" means using objects where objects fit — not making everything a class.

## Modern syntax used throughout

| Syntax | Meaning | Where |
|---|---|---|
| `?.` | optional chaining — `a?.b` is `undefined` if `a` is nullish | `row.month_day?.trim()` |
| `??` | nullish coalescing — falls back only on `null`/`undefined`, unlike `\|\|` | `region ?? 'Unknown'` |
| `??=` | assign only if nullish | `sharedDb ??= new PostgresDatabase(...)` |
| `...` | spread / rest | `{ ...state, ...patch }` |
| `{ a, b }` | destructuring | `constructor({ id, code, name })` |
| `` `${}` `` | template literals | error messages |

`??` versus `||` matters: `0 || 5` is `5`, but `0 ?? 5` is `0`. With event counts
and years, that difference is a real bug.

## Do it yourself

1. **Feel the missing compiler.** In `TimelineService`, swap the parameters of
   `getEventsForYearAndCountry(year, countryCode)` without touching the
   controller. Run `npm test`. Note that the failure comes from a *test*, not a
   compiler — and that without the test it would have reached a user.

2. **Add an error type.** Create `ConflictError extends DomainError` with
   `code = 'CONFLICT'` and `httpStatus = 409`. Throw it from a service method.
   Confirm the API returns 409 **without editing `errorHandler.js`**.

3. **Break Liskov deliberately.** Delete `.sort(...)` from
   `InMemoryCountryRepository.findAll()`. Run
   `node --test server/tests/repositoryContract.test.js`. Watch it get caught.

4. **Meet the `this` bug.** In `TimelineController.routes()`, change
   `asyncRoute(this.#getTimeline)` to a normal method reference and watch it
   break. Understand *why* before you undo it.

## Check yourself

- Why is `Year` immutable, and what would break if it were not?
- Why does `Country` compare by `id` while `Year` compares by `value`?
- Why does `CountryRepository`'s constructor check `new.target`?
- What is the difference between `??` and `||`?

<details>
<summary>Answers</summary>

- **Immutability**: a `Year` can be shared and cached safely. If a holder could
  mutate it, a `Year` passed into a function could come back different — and
  `#private` with no setter means that is not expressible.
- **`Country` vs `Year`**: `Country` is an *entity* — identity survives its
  attributes changing (fix a typo in the name, still the same country). `Year` is
  a *value object* — 1969 is 1969; there is no "which 1969". Chapter 3 is
  entirely about this.
- **`new.target`**: it is how you express "abstract" in a language without the
  keyword. `new CountryRepository()` is always a mistake, and failing loudly at
  that moment beats a confusing error later.
- **`??` vs `||`**: `||` falls back on any falsy value (`0`, `''`, `false`); `??`
  only on `null`/`undefined`. For `eventCount ?? 0`, using `||` would be
  harmless, but for a *legitimate* `0` you want to keep, `||` silently replaces
  it.

</details>

→ Next: [03 — The domain model](03-domain-model.md)
