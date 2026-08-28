import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import { Store, Async } from '../js/Store.js';

// The Store touches no DOM, so this file needs no jsdom at all — one more sign
// that keeping state separate from rendering was the right split.
describe('Store', () => {
  it('exposes its initial state', () => {
    assert.deepEqual(new Store({ a: 1 }).getState(), { a: 1 });
  });

  it('merges a patch rather than replacing everything', () => {
    const store = new Store({ a: 1, b: 2 });
    store.setState({ b: 3 });
    assert.deepEqual(store.getState(), { a: 1, b: 3 });
  });

  it('replaces the state object instead of mutating it', () => {
    const store = new Store({ a: 1 });
    const before = store.getState();
    store.setState({ a: 2 });
    // The old object is untouched, so `previous !== next` is a reliable way to
    // detect a change — which is what the views rely on.
    assert.equal(before.a, 1);
    assert.notEqual(before, store.getState());
  });

  it('freezes state so a view cannot mutate it by accident', () => {
    const store = new Store({ a: 1 });
    assert.throws(() => {
      store.getState().a = 99;
    }, TypeError);
  });

  it('notifies subscribers with the new and previous state', () => {
    const store = new Store({ a: 1 });
    const calls = [];
    store.subscribe((next, previous) => calls.push([next.a, previous.a]));
    store.setState({ a: 2 });
    assert.deepEqual(calls, [[2, 1]]);
  });

  it('does not notify when nothing actually changed', () => {
    const store = new Store({ a: 1 });
    let calls = 0;
    store.subscribe(() => (calls += 1));
    store.setState({ a: 1 });
    assert.equal(calls, 0);
  });

  it('returns an unsubscribe function', () => {
    const store = new Store({ a: 1 });
    let calls = 0;
    const unsubscribe = store.subscribe(() => (calls += 1));
    store.setState({ a: 2 });
    unsubscribe();
    store.setState({ a: 3 });
    assert.equal(calls, 1);
  });

  it('supports several independent subscribers', () => {
    const store = new Store({ a: 1 });
    let one = 0;
    let two = 0;
    store.subscribe(() => (one += 1));
    store.subscribe(() => (two += 1));
    store.setState({ a: 2 });
    assert.equal(one, 1);
    assert.equal(two, 1);
  });

  it('survives a subscriber unsubscribing during notification', () => {
    const store = new Store({ a: 1 });
    let secondRan = false;
    const off = store.subscribe(() => off());
    store.subscribe(() => (secondRan = true));
    assert.doesNotThrow(() => store.setState({ a: 2 }));
    assert.ok(secondRan);
  });

  it('fails loudly instead of looping if a subscriber calls setState', () => {
    const store = new Store({ a: 1 });
    store.subscribe(() => store.setState({ a: Math.random() }));
    assert.throws(() => store.setState({ a: 2 }), /while notifying/);
  });
});

describe('Async state helpers', () => {
  it('models the four states', () => {
    assert.equal(Async.idle().status, 'idle');
    assert.equal(Async.loading().status, 'loading');
    assert.equal(Async.success([1]).status, 'success');
    assert.equal(Async.error(new Error('x')).status, 'error');
  });

  it('carries data only in success and error only in error', () => {
    assert.deepEqual(Async.success(['a']).data, ['a']);
    assert.equal(Async.success(['a']).error, null);
    assert.equal(Async.error(new Error('boom')).data, null);
    assert.match(Async.error(new Error('boom')).error.message, /boom/);
  });
});
