import { describe, it, before } from 'node:test';
import assert from 'node:assert/strict';
import { setupDom } from './domEnvironment.js';

before(() => setupDom());

const { el, append, clear, render } = await import('../js/dom.js');

describe('el', () => {
  it('creates an element with text', () => {
    const node = el('p', {}, 'hello');
    assert.equal(node.tagName, 'P');
    assert.equal(node.textContent, 'hello');
  });

  it('sets class, properties and dataset', () => {
    const node = el('a', { class: 'x y', href: 'https://example.org/', dataset: { year: '1969' } });
    assert.equal(node.className, 'x y');
    assert.equal(node.getAttribute('href'), 'https://example.org/');
    assert.equal(node.dataset.year, '1969');
  });

  it('sets ARIA and role as real attributes', () => {
    const node = el('button', { role: 'radio', 'aria-checked': 'true', tabindex: -1 });
    assert.equal(node.getAttribute('role'), 'radio');
    assert.equal(node.getAttribute('aria-checked'), 'true');
    assert.equal(node.getAttribute('tabindex'), '-1');
  });

  it('attaches event listeners', () => {
    let clicks = 0;
    const node = el('button', { on: { click: () => (clicks += 1) } });
    node.click();
    node.click();
    assert.equal(clicks, 2);
  });

  it('nests children and flattens arrays', () => {
    const node = el('ul', {}, [el('li', {}, 'a'), [el('li', {}, 'b'), el('li', {}, 'c')]]);
    assert.equal(node.querySelectorAll('li').length, 3);
  });

  it('skips null and false children, so `cond && el(...)` works', () => {
    const node = el('div', {}, [null, false, undefined, el('span', {}, 'kept')]);
    assert.equal(node.childNodes.length, 1);
  });

  it('skips null and false props', () => {
    const node = el('button', { 'aria-current': null, disabled: false });
    assert.equal(node.hasAttribute('aria-current'), false);
    assert.equal(node.disabled, false);
  });
});

/**
 * THE SECURITY TEST.
 *
 * Event titles and summaries come from a database. If the DOM helper ever
 * started parsing them as HTML, a malicious entry would become executable
 * markup — cross-site scripting. This test pins that it does not, and it is the
 * kind of test worth keeping forever, because the unsafe version LOOKS FINE
 * until the day someone exploits it.
 */
describe('el is XSS-safe', () => {
  const attack = '<img src=x onerror="globalThis.__pwned = true">';

  it('renders markup in text as literal characters, not as elements', () => {
    const node = el('h3', {}, attack);
    assert.equal(node.querySelector('img'), null, 'must not create a real <img>');
    assert.equal(node.textContent, attack);
    assert.equal(globalThis.__pwned, undefined);
  });

  it('is safe for nested children too', () => {
    const node = el('div', {}, [el('p', {}, attack)]);
    assert.equal(node.querySelector('img'), null);
  });

  it('demonstrates the unsafe alternative for contrast', () => {
    // This is what `innerHTML = ...` would have done: real markup.
    const unsafe = document.createElement('h3');
    unsafe.innerHTML = attack;
    assert.notEqual(unsafe.querySelector('img'), null, 'innerHTML parses it as an element');
  });
});

describe('clear and render', () => {
  it('clear empties a node', () => {
    const node = el('div', {}, [el('span'), el('span')]);
    clear(node);
    assert.equal(node.childNodes.length, 0);
  });

  it('render replaces the contents', () => {
    const node = el('div', {}, el('span', {}, 'old'));
    render(node, el('p', {}, 'new'));
    assert.equal(node.textContent, 'new');
    assert.equal(node.querySelectorAll('span').length, 0);
  });

  it('append returns the parent so calls can be chained', () => {
    const node = el('div');
    assert.equal(append(node, 'x'), node);
  });
});
