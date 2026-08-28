# 03 — The domain model

> **Concept** → The code → Do it yourself → Check yourself

The domain layer is `server/src/domain/`. Four files, no imports from anywhere
else in the project, and every rule the product cares about.

## Entities and value objects

This is the distinction that makes domain modelling click.

**A value object has no identity.** It *is* its value. Two `Year`s with the value
1969 are interchangeable — asking "which 1969?" is meaningless. They compare by
value, and they are immutable.

**An entity has identity that outlives its attributes.** Country #12 is South
Africa. Fix a typo in the name, change its region — still country #12.

```js
// Value object: compare the contents.
equals(other) { return other instanceof Year && other.value === this.#value; }

// Entity: compare the identity.
equals(other) { return other instanceof Country && other.id === this.id; }
```

Getting this backwards is a classic bug: compare two entities field-by-field and
a renamed country suddenly looks like a *different* country, silently
duplicating rows.

**How to tell them apart:** ask *"if two of these have identical fields, are they
the same thing?"* Yes → value object. No → entity.

| | `Year` | `Country` | `HistoricalEvent` |
|---|---|---|---|
| Kind | value object | entity | entity |
| Compares by | value | id | id |
| Immutable | yes (`#private`) | yes (`freeze`) | yes (`freeze`) |
| Has an `id` | no | yes | yes |

## Making bad states unrepresentable

The `Year` constructor rejects non-numbers, non-integers, anything before 100 AD
and anything in the future. So this is *impossible* anywhere in the system:

```js
new Year(3050);      // throws ValidationError
new Year('1969');    // throws — a string is not a number
```

There is no code path that produces a `Year` of 3050. Not from the API, not from
the seed script, not from a test.

Compare that with passing plain numbers around and checking
`if (year > currentYear)` in — how many places? You will never be sure you got
them all. **Push validation to the moment of construction, then let the type
carry the guarantee.**

`HistoricalEvent` takes this one step further:

```js
if (!(year instanceof Year)) {
  throw new ValidationError('Event year must be a Year instance.', 'year');
}
```

Accepting a raw number here would let an unvalidated year in through the back
door, defeating the point of having a `Year` class at all. In a typed language
the compiler would enforce this; in JavaScript you write the check.

## The editorial rule as code

The clearest example in this project of a *product* rule becoming an
*engineering* rule:

```js
if (typeof sourceUrl !== 'string' || !/^https?:\/\/\S+$/.test(sourceUrl)) {
  throw new ValidationError(`Event must cite a source URL, ...`, 'sourceUrl');
}
```

"This is a factual, educational history site, so every claim cites a source" is
an *editorial* promise. Here it is a fact about the program. An uncited event
cannot be constructed — so it cannot be seeded, cannot be returned by the API,
cannot be rendered, and cannot appear in a test fixture.

The test that pins it (`tests/domain.test.js`):

```js
it('refuses to exist without a citable source', () => {
  for (const sourceUrl of ['', 'not-a-url', 'ftp://example.com/x', 'wikipedia.org/wiki/X', null]) {
    assert.throws(() => anEvent({ sourceUrl }), ValidationError);
  }
});
```

**Find the equivalent rule in whatever you build next.** Every product has one —
the promise that, if broken, means the product is not the thing you said it was.
Put it in a constructor.

## Modelling imprecise history

A design decision worth dwelling on, because a naive model gets it wrong:

```js
this.year = year;              // a Year — always known
this.monthDay = monthDay;      // "04-27", or null
```

Why not a `Date`? Because **`Date` cannot represent "sometime in 1347"**, and
most of history is like that. A `Date` would force you to invent a day — and an
invented 1 January is a lie the system can never distinguish from a real one.

Two fields, one required and one optional, model *the actual state of historical
knowledge*. The ordering logic follows (`compareByDate`): dated events first, in
date order; undated events after, alphabetically. And the UI shows "27 April
1994" or just "1994" — never a fabricated date.

> The general lesson: **model what you actually know, including the gaps.** When
> your model cannot express uncertainty, your code will invent certainty.

## An enum without `enum`

JavaScript has no `enum`, so:

```js
export const EventCategory = Object.freeze({
  Politics: 'politics',
  Conflict: 'conflict',
  ...
});
export const EVENT_CATEGORIES = Object.freeze(Object.values(EventCategory));
```

`Object.freeze` means a typo like `EventCategory.Politcs` evaluates to
`undefined` — which then fails the category check — rather than silently adding
a new key to the object.

## The error hierarchy

```
Error (built-in)
└── DomainError (abstract) ── code, httpStatus
    ├── ValidationError ──── VALIDATION_ERROR, 400
    └── NotFoundError ────── NOT_FOUND, 404
```

`DomainError` refuses direct construction via `new.target`, so every error is a
specific kind. Each subclass carries its own `httpStatus`, which lets the API
layer map errors polymorphically instead of with a growing `if/else`.

What this buys: **`throw` in the domain, and the right HTTP status comes out the
other end.** The service never mentions 404. The controller never mentions
`NotFoundError`.

## Serialisation: `toJSON` as a boundary

```js
toJSON() {
  return { id: this.id, code: this.code, name: this.name, region: this.region };
}
```

Explicit, not automatic. `JSON.stringify` would dump every public property — so
the day someone adds `internalNotes` to `Country`, it silently ships to every
browser. An explicit `toJSON` is a **deliberate boundary between your domain and
the outside world**, and the integration test asserts on it:

```js
assert.deepEqual(Object.keys(body.countries[0]).sort(), ['code', 'id', 'name', 'region']);
```

That test fails if a field is ever added or removed — exactly what you want from
a public contract.

(`Year.toJSON()` returns a bare number, which is why `event.year` arrives in the
browser as `1969` rather than as an object.)

## Do it yourself

1. **Add a category.** Add `Migration: 'migration'` to `EventCategory`. Then find
   everything that needs to change: the SQL `CHECK` constraint needs a new
   migration, and the CSS wants a colour. Notice that *nothing* found this for
   you — in a typed language the compiler would have. Write a test that asserts
   every `EVENT_CATEGORIES` value has a CSS rule, and you have replaced the
   compiler with a test.

2. **Add an invariant.** Make `HistoricalEvent` reject a summary shorter than 20
   characters ("a one-word summary teaches nothing"). Write the test first, watch
   it fail, then make it pass. Then run the full suite — does the seed data
   survive? `seedData.test.js` will tell you immediately.

3. **Model a new value object.** Events often have a *place* — "Amritsar", "Kitty
   Hawk". Design a `Place`: what makes it valid? Is it a value object or an
   entity? (There is a case for either, and the answer depends on whether two
   places with the same name are the same place.)

## Check yourself

- Why does `Year` have a static `fromString` as well as a constructor?
- What would break if `HistoricalEvent` were mutable?
- Why does the database *also* enforce the source-URL rule?

<details>
<summary>Answers</summary>

- **`fromString`**: the HTTP boundary hands you strings, and string parsing is a
  separate concern from year validation. A named factory documents intent and
  keeps the constructor focused. `Year.fromString('abc')` fails with *"not a
  valid year"*; `new Year(NaN)` fails with *"must be a whole number"* — two
  different mistakes, two different messages.
- **Mutability**: an event handed to two views could be changed by one and
  observed changed by the other — a bug with no stack trace pointing at the
  culprit. `Object.freeze` turns silent corruption into an immediate throw.
- **Belt and braces**: application code can be bypassed. A colleague running
  `psql` at 2am, a bulk import, a future service in another language — none of
  them go through your constructor. The database constraint is the only defence
  that is truly unavoidable.

</details>

→ Next: [04 — PostgreSQL](04-postgres.md)
