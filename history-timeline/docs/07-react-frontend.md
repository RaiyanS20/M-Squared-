# 07 — The React frontend

> **Concept** → The code → Do it yourself → Check yourself

## The mental model

React is one idea: **UI is a function of state.** You do not update the DOM; you
change state and describe what the UI looks like for that state. React works out
the DOM changes.

```
state ──▶ render ──▶ what the user sees
  ▲                        │
  └──────── events ────────┘
```

Everything else — hooks, effects, keys — is machinery serving that idea.

## State design: the decision that matters most

`App.tsx` stores exactly **two** things:

```ts
const [selectedYear, setSelectedYear] = useState<number | null>(null);
const [selectedCountryCode, setSelectedCountryCode] = useState<string | null>(null);
```

Everything else on screen is **derived**:

```ts
const timeline  = useTimeline();                                  // from nothing
const countries = useCountriesForYear(selectedYear);              // from year
const events    = useEvents(selectedYear, selectedCountryCode);   // from both
```

There is no `countries` state, no `events` state, no `isLoading` flag. So there
is **no way for the panels to disagree with each other** — the events shown are
always the events for the year and country in the heading, because both come
from the same two values.

> **The rule: if you can compute it, do not store it.** Almost every "the UI is
> showing stale data" bug is a stored copy of something derivable.

## Why `null` is a real state

`selectedYear: number | null` — `null` means "nothing chosen yet", which is
genuinely different from any year. The types force every consumer to handle it,
which is why the second panel says *"Select a year to begin"* instead of
rendering an empty box.

## The three states of remote data — and the fourth

`useAsyncResource` models them as a discriminated union:

```ts
export type AsyncState<T> =
  | { status: 'idle';    data: null; error: null }
  | { status: 'loading'; data: null; error: null }
  | { status: 'success'; data: T;    error: null }
  | { status: 'error';   data: null; error: ApiError };
```

Compare the common alternative — `const [data, setData] = useState(); const
[loading, setLoading] = useState(); const [error, setError] = useState();` —
which permits `loading: true` *and* `error` set *and* stale `data`, all at once.
Eight combinations, of which four are nonsense.

With a union, **the impossible states cannot be written down**, and TypeScript
narrows `data` to non-null inside the `success` branch:

```ts
{events.status === 'success' && <EventList events={events.data.events} ... />}
//                                                     ^^^^ known non-null here
```

### The `idle` state exists because of a real bug

`idle` was added while writing this project, after a test failed. The App test
*"clears the chosen country when the year changes"* caught it:

Changing the year clears `selectedCountryCode`, which disables the `useEvents`
resource. The hook originally just *stopped*, leaving its last `success` state
in place — so **the previous country's events stayed on screen under the new
year's heading.** Plausible-looking, completely wrong data.

The fix is in the hook, not in the component:

```ts
if (!enabled) {
  dispatch({ type: 'idle' });   // a disabled resource forgets what it loaded
  return;
}
```

Two lessons. First: **a disabled resource must forget its data.** Second, and
bigger: that bug was invisible by eye — the screen looked fine — and a test found
it. This is what tests are actually for.

## Race conditions, and the AbortController

Click 1914, then quickly click 1969. Two requests are in flight. If 1914's
response happens to arrive *second*, it overwrites 1969's — and the screen now
shows 1914's countries under a 1969 heading.

This is invisible on a fast local network and very visible on a phone.

```ts
useEffect(() => {
  const controller = new AbortController();
  dispatch({ type: 'loading' });

  load(controller.signal)
    .then((data) => { if (!controller.signal.aborted) dispatch({ type: 'success', data }); })
    .catch((error) => { if (controller.signal.aborted) return; ... });

  return () => controller.abort();   // cleanup: runs before the next effect
}, [...deps, enabled]);
```

The cleanup function runs **before the effect re-runs** and on unmount. So the
old request is cancelled the moment a new one starts. The bug is not made
*unlikely*; it is made **impossible**.

`ApiClient` cooperates by re-throwing `AbortError` untouched rather than wrapping
it — an abort is a normal part of the lifecycle, not a failure to report:

```ts
if (cause instanceof DOMException && cause.name === 'AbortError') throw cause;
```

**Every `useEffect` that starts something must clean it up.** Requests,
subscriptions, timers, listeners. React's StrictMode double-invokes effects in
development specifically to make missing cleanup obvious.

## Controlled components

`Timeline` holds no state:

```tsx
<Timeline ticks={...} selectedYear={selectedYear} onSelectYear={handleSelectYear} />
```

It receives what to show and reports what was clicked. All state lives in `App`.
This is *"lifting state up"*, and it is why the timeline, the country list and
the heading can never disagree.

A component holding its own copy of `selectedYear` would be a second source of
truth — and two sources of truth eventually differ.

## Dependent state must be reset

```ts
const handleSelectYear = useCallback((year: number) => {
  setSelectedYear(year);
  setSelectedCountryCode(null);   // ← the important line
}, []);
```

"France in 1969" is valid; after switching to 1848, France may have nothing
recorded. **When a parent selection changes, clear what depended on it.** This
one line prevents a whole family of stale-UI bugs — and combined with the `idle`
fix above, the events panel correctly returns to its prompt.

## Accessibility, done properly

The timeline is not a row of clickable `div`s. It is a real ARIA radiogroup:

```tsx
<div role="radiogroup" aria-label="Choose a year" onKeyDown={handleKeyDown}>
  <button role="radio" aria-checked={isSelected} tabIndex={isSelected ? 0 : -1}
          aria-label={`${tick.year}, ${tick.eventCount} events`}>
```

- **`<button>`, not `<div>`** — keyboard operable and announced as a control.
- **Roving tabindex** — one year is in the tab order; arrows move between them.
  Tabbing through 127 buttons is not navigation, it is a punishment.
- **`aria-label` with the count** — a screen reader says *"1994, 1 event"*. The
  bar height is a visual affordance; the label carries the same information
  non-visually.
- **Arrow / PageUp / PageDown / Home / End**, clamped at both ends.

The event cards' links carry `rel="noopener noreferrer"`: `noopener` stops the
opened page reaching back through `window.opener`, `noreferrer` withholds the
referring URL.

None of this is optional polish. An educational site that a keyboard or screen
reader user cannot operate has failed at being educational.

## Classes on the frontend: where they still belong

Components are functions. The **network layer** is classes:

```ts
export abstract class ApiClient {
  protected async get<T>(path: string, signal?: AbortSignal): Promise<T> { ... }
}
export class TimelineApi extends ApiClient {
  getTimeline(startYear?, endYear?, signal?) { return this.get(...); }
}
```

`ApiClient` owns *every* HTTP concern — URL building, headers, non-2xx handling,
malformed JSON, network failure, aborts — so each `TimelineApi` method is one
typed line. Adding authentication or retries is a change to one file.

**No component ever calls `fetch`.** Same principle as services never writing SQL.

## Deriving with `useMemo`

```ts
export function useCountriesByRegion<T extends { region: string; name: string }>(
  countries: readonly T[] | undefined,
): Array<[string, T[]]>
```

Grouping runs only when `countries` changes, not on every render.

**Do not reach for `useMemo` by default.** It has a cost — the comparison and the
retained reference. Use it for genuinely expensive work or to stabilise a
reference that something else depends on. Grouping a dozen countries is cheap;
this one is as much about intent as speed.

## CSS without a framework

Plain CSS, custom properties, BEM-ish naming (`block__element--modifier`). No
Tailwind, no CSS-in-JS — because the goal is to learn what those tools *do for
you*. Custom properties, flexbox and grid cover this entire app.

Two details worth stealing:

**Dark mode is a token swap.** Colours are declared once as custom properties and
redefined inside `@media (prefers-color-scheme: dark)`. No component knows a
theme exists.

**A specificity bug, and its fix.** The selected country was rendering
light-on-light while hovered:

```css
.countries__item:hover          { }  /* specificity 0,2,0 — class + pseudo-class */
.countries__item--selected      { }  /* specificity 0,1,0 — loses */
```

The fix is to exclude the selected row from the hover rule rather than escalate:

```css
.countries__item:not(.countries__item--selected):hover { background: var(--surface-alt); }
```

`!important` would have "worked" and left a landmine. **Understand the cascade;
do not fight it.**

Also honoured: `prefers-reduced-motion` (animation makes some people ill) and a
visible `:focus-visible` ring (it is how keyboard users navigate).

## Do it yourself

1. **See the race condition.** In `ApiClient.get`, add a random delay before
   `fetch`:
   ```ts
   await new Promise((r) => setTimeout(r, Math.random() * 2000));
   ```
   Click several years quickly. Then delete the `return () => controller.abort()`
   line from `useAsyncResource` and try again. Watch the panels disagree. Restore
   both.

2. **Reproduce the stale-data bug.** Remove the `dispatch({ type: 'idle' })`
   branch from `useAsyncResource`. Run `npx vitest run App`. Read the failure —
   it is the exact bug described above.

3. **Test the keyboard.** Load the app, press Tab to the timeline, then arrows,
   PageUp/PageDown, Home/End. Now do it with your OS screen reader on.

4. **Add a feature.** Filter events by category. Where does the state live —
   `App`, or `EventList`? (Consider: should it survive changing country?)

5. **Deep-link it.** Put `selectedYear` and `selectedCountryCode` in the URL
   (`/1994/ZAF`) so a link can be shared. This is the natural moment to introduce
   a router, and it will make you notice how much the two-state design helped.

## Check yourself

- Why does `useAsyncResource` take `deps` explicitly instead of depending on
  `load`?
- Why is the timeline a controlled component?
- Why is `status` a union rather than three booleans?

<details>
<summary>Answers</summary>

- **Explicit `deps`**: callers pass an inline arrow function, which is a new
  reference every render. Including `load` in the dependency array would re-run
  the effect forever. `deps` is the honest list of what the request actually
  depends on. (The alternative is wrapping every loader in `useCallback` — more
  ceremony at every call site for the same result.)
- **Controlled**: one source of truth. Local state in the timeline could drift
  from `App`'s, and then the heading and the highlighted year disagree.
- **Union over booleans**: three booleans allow eight combinations, half of them
  nonsense (`loading && error`). A union permits only the four real states, and
  lets TypeScript narrow `data` to non-null in the success branch — so
  `data.events` needs no `?.`.

</details>

→ Next: [08 — Testing](08-testing.md)
