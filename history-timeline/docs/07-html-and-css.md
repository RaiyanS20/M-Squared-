# 07 — HTML and CSS

> **Concept** → The code → Do it yourself → Check yourself

Two files: `client/index.html` and `client/css/styles.css`. No framework, no
preprocessor, no build step. Everything here is the platform.

## HTML is structure, not appearance

`index.html` describes **what the page is**, never how it looks:

```html
<main class="app__main">
  <section class="panel" aria-labelledby="timeline-heading">
    <h2 id="timeline-heading" class="panel__heading">1. Choose a year</h2>
    <div id="timeline"></div>
  </section>
</main>
```

`<section>`, `<main>`, `<header>`, `<h2>`, `<nav>`, `<article>`, `<time>`,
`<ol>` — these are **semantic** elements. They carry meaning that `<div>` does
not.

**Why it matters, concretely:**

- A screen reader user can jump between headings and landmarks. With
  `<div class="heading">` there is nothing to jump to.
- Search engines and preview cards read the structure.
- Browser reader modes work.
- Your CSS gets simpler, because the structure already says what things are.

The rule: **choose the element that describes the content, then style it.** Reach
for `<div>` only when nothing else fits — it means "a box with no meaning", and
that is occasionally the honest answer.

### Structure in HTML, data in JavaScript

The skeleton is written by hand, once. JavaScript fills in three regions:

```html
<div id="timeline"></div>
<div id="countries"></div>
<div id="events"></div>
```

The page is meaningful before a single line of JavaScript runs, and the
JavaScript stays focused on data rather than re-declaring the page shape on every
render.

### `<script type="module">`

```html
<script type="module" src="/js/app.js"></script>
```

Four things at once: `import`/`export` work natively, the script is deferred (so
the DOM exists when it runs), it runs in strict mode, and its top-level variables
do not leak onto `window`. **This one attribute is why the project needs no
bundler.**

### Accessibility attributes you will actually use

```html
<section aria-labelledby="timeline-heading">   <!-- names the region -->
<span class="visually-hidden">opens in a new tab</span>
```

And from `TimelineView.js`:

```html
<div role="radiogroup" aria-label="Choose a year">
  <button role="radio" aria-checked="true" tabindex="0" aria-label="1994, 1 event">
```

**Use a real `<button>`.** A `<div onclick>` is not focusable, not reachable by
keyboard, and not announced as a control. Every interactive thing in this project
is a `<button>` or an `<a>`, and that decision alone does most of the
accessibility work.

The distinction: **`<a>` navigates, `<button>` acts.** The source links are `<a>`
because they go somewhere; everything else is a `<button>`.

## CSS: the cascade, specificity, and one real bug

CSS resolves conflicts by **specificity**, counted as (ids, classes, elements):

| Selector | Specificity | |
|---|---|---|
| `p` | 0,0,1 | element |
| `.panel` | 0,1,0 | class |
| `.item:hover` | 0,2,0 | class + **pseudo-class counts as a class** |
| `#timeline` | 1,0,0 | id |

Higher wins. Ties go to whichever comes **last** in the file.

**This project hit exactly that bug.** The selected country rendered
light-on-light while the mouse rested on it:

```css
.countries__item:hover     { background: var(--surface-alt); }  /* 0,2,0 — wins */
.countries__item--selected { background: var(--accent); }       /* 0,1,0 — loses */
```

The fix is to *exclude* rather than escalate:

```css
.countries__item:not(.countries__item--selected):hover { background: var(--surface-alt); }
```

`!important` would also have "worked" and left a landmine for the next person.
**Understand the cascade; do not fight it.**

This is also why the class names look like `block__element--modifier` (BEM). It
is a naming convention, not a technology, and its entire purpose is to keep every
selector at one class so specificity stays flat and predictable.

## Custom properties (CSS variables)

```css
:root {
  --accent: #8a5a2b;
  --ink: #23201b;
}
.event__source { color: var(--accent); }
```

Declared once, used everywhere — changing the accent colour of the whole app is
one line.

Unlike a preprocessor variable, these are **live in the browser**: they cascade,
they can be read and changed from JavaScript, and they can be redefined inside a
media query. Which gives us dark mode almost free:

```css
@media (prefers-color-scheme: dark) {
  :root { --bg: #17150f; --ink: #ece7dc; --accent: #d8a05e; }
}
```

**No component knows a theme exists.** They all use `var(--ink)`, and the tokens
change underneath. That is the whole implementation of dark mode in this project
— and it respects the reader's operating-system setting rather than imposing a
choice.

## Layout: flexbox and grid

Two systems, and the choice between them is genuinely simple:

- **Flexbox** — content laid out along *one* axis. Use it for a row or a column.
- **Grid** — content laid out in *two* dimensions at once, or when you want to
  define the tracks explicitly.

```css
/* Flexbox: the timeline is a single row of bars. */
.timeline__rail { display: flex; align-items: flex-end; gap: 2px; overflow-x: auto; }

/* Grid: two columns side by side. */
.app__columns { display: grid; grid-template-columns: minmax(240px, 1fr) 2fr; }
```

`gap` works in both and has replaced margin hacks for spacing between items.

### Mobile first

```css
.app__columns { grid-template-columns: 1fr; }          /* the default: one column */

@media (min-width: 820px) {                            /* the enhancement */
  .app__columns { grid-template-columns: minmax(240px, 1fr) 2fr; }
}
```

The single-column layout is the **default**, and the wide layout is added inside
a `min-width` query. Written the other way round you end up overriding desktop
rules on the device with the least power and the smallest screen.

`clamp(1rem, 3vw, 2.5rem)` gives fluid padding — never below 1rem, never above
2.5rem — with no media query at all.

## Details that separate a finished page from a demo

**Never remove the focus ring.**

```css
:focus-visible { outline: 3px solid var(--focus); outline-offset: 2px; }
```

`:focus-visible` shows it for keyboard users but not on mouse clicks, which is
why you can style it boldly. `outline: none` with no replacement makes a site
unusable by keyboard.

**Honour reduced motion.** Animation makes some people physically ill, and their
OS setting says so:

```css
@media (prefers-reduced-motion: reduce) {
  *, *::before, *::after {
    animation-duration: 0.01ms !important;
    transition-duration: 0.01ms !important;
  }
}
```

(One of the few legitimate uses of `!important`: overriding everything is the
literal intent.)

**Colour must never be the only signal.** Event categories are colour-coded *and*
labelled in text, so the meaning survives colour blindness and monochrome print.

**Buttons do not inherit your font.** `font: inherit` on every button, or they
render in the browser's default control font and look out of place.

**Tabular numerals** — `font-variant-numeric: tabular-nums` — stop year labels
jittering as digits change width.

## Do it yourself

1. **Retheme the app in one line.** Change `--accent` in `:root` and reload.
   Then change it only inside the dark-mode block.

2. **Reproduce the specificity bug.** Change the hover rule back to
   `.countries__item:hover` and watch the selected country become unreadable
   while hovered. Fix it three ways — `:not()`, reordering, `!important` — and
   decide which you would want to inherit.

3. **Break the semantics.** Change the `<button>` elements in `TimelineView.js`
   to `<div>`. The mouse still works. Now try Tab and the arrow keys. That gap is
   what semantic HTML gives you for free.

4. **Test it at 320px wide** in your browser's device toolbar. Then zoom the page
   to 200% — text should reflow, not clip.

5. **Read the page with a screen reader** (VoiceOver on macOS: ⌘F5; Narrator on
   Windows: Ctrl+Win+Enter). Navigate by headings. This is the single most
   clarifying hour you can spend on front-end work.

## Check yourself

- Why `<section aria-labelledby="...">` rather than just a `<div>`?
- Why is `--accent` defined on `:root` rather than on each component?
- Why does the mobile layout come first?

<details>
<summary>Answers</summary>

- **`<section>` + `aria-labelledby`**: it creates a named landmark, so assistive
  technology can list and jump between the page's regions. A `<div>` creates
  nothing to navigate to.
- **Tokens on `:root`**: `:root` is the document root, so the variable cascades
  everywhere and can be redefined once (in a media query, or by JavaScript) to
  restyle the whole app. Defined per component you would have forty places to
  change and no way to swap a theme.
- **Mobile first**: `min-width` queries mean the simplest layout is the default
  and complexity is added as space allows. The reverse forces the least capable
  device to parse and override the most complex rules.

</details>

→ Next: [08 — JavaScript in the browser](08-javascript-in-the-browser.md)
