import { el, render } from '../dom.js';
import { loadingMessage, errorMessage } from './statusMessages.js';

/**
 * MVP features 1 and 2: the scale, and choosing a year on it.
 *
 * ---------------------------------------------------------------------------
 * THE MOST IMPORTANT FILE IN THE COURSE FOR UNDERSTANDING FRAMEWORKS.
 * ---------------------------------------------------------------------------
 *
 * A naive version would rebuild all 126 buttons every time the state changed.
 * It would work, and it would be wrong in ways you can feel:
 *
 *   * the focused button is destroyed mid-keypress, so keyboard navigation dies
 *   * the horizontal scroll position jumps back to the start
 *   * 126 elements are discarded and recreated to change one CSS class
 *
 * So this view separates two jobs:
 *
 *   #renderRail()      builds the buttons ONCE, when the data arrives, and
 *                      remembers each one in a Map keyed by year.
 *   #applySelection()  changes only the two buttons that actually differ.
 *
 * That bookkeeping — "what changed, and what is the minimum DOM edit?" — is
 * precisely what React's virtual DOM does for you automatically. Doing it by
 * hand once tells you what you are buying when you later choose a framework,
 * and why the answer is not always "yes".
 */
export class TimelineView {
  #root;
  #onSelectYear;
  #buttons = new Map();
  #rail = null;
  #renderedTicks = null;
  #selectedYear = null;

  constructor(root, { onSelectYear }) {
    this.#root = root;
    this.#onSelectYear = onSelectYear;
  }

  /** Called on every state change. Decides how little work it can get away with. */
  update(state) {
    const { timeline, selectedYear } = state;

    if (timeline.status === 'loading' || timeline.status === 'idle') {
      this.#reset();
      render(this.#root, loadingMessage('Loading the timeline…'));
      return;
    }

    if (timeline.status === 'error') {
      this.#reset();
      render(this.#root, errorMessage(timeline.error));
      return;
    }

    // Rebuild the rail only when the DATA changed — not when the selection did.
    if (timeline.data.ticks !== this.#renderedTicks) {
      this.#renderRail(timeline.data.ticks);
    }

    if (selectedYear !== this.#selectedYear) {
      this.#applySelection(selectedYear);
    }
  }

  #reset() {
    this.#buttons.clear();
    this.#rail = null;
    this.#renderedTicks = null;
    this.#selectedYear = null;
  }

  #renderRail(ticks) {
    this.#buttons.clear();
    this.#selectedYear = null;

    // Bar heights are relative to the busiest year in view. `Math.max(1, ...)`
    // avoids dividing by zero when nothing is recorded at all.
    const busiest = Math.max(1, ...ticks.map((t) => t.eventCount));

    const buttons = ticks.map((tick) => {
      const isDecade = tick.year % 10 === 0;
      const noun = tick.eventCount === 1 ? 'event' : 'events';

      const button = el(
        'button',
        {
          type: 'button',
          // A real ARIA radiogroup. A row of clickable <div>s is invisible to
          // assistive technology and unreachable from a keyboard.
          role: 'radio',
          'aria-checked': 'false',
          tabindex: -1,
          // The visible bar carries this information to sighted users; the label
          // carries the same information to everyone else.
          'aria-label': `${tick.year}, ${tick.eventCount} ${noun}`,
          title: `${tick.year} — ${tick.eventCount} recorded ${noun}`,
          class: `timeline__year${isDecade ? ' timeline__year--decade' : ''}${
            tick.eventCount === 0 ? ' timeline__year--empty' : ''
          }`,
          dataset: { year: String(tick.year) },
          on: { click: () => this.#onSelectYear(tick.year) },
        },
        [
          el('span', {
            class: 'timeline__bar',
            'aria-hidden': 'true',
            style: { height: `${(tick.eventCount / busiest) * 100}%` },
          }),
          el('span', { class: 'timeline__label', 'aria-hidden': 'true' }, isDecade ? String(tick.year) : ''),
        ],
      );

      this.#buttons.set(tick.year, button);
      return button;
    });

    this.#rail = el(
      'div',
      {
        class: 'timeline__rail',
        role: 'radiogroup',
        'aria-label': 'Choose a year',
        // ONE listener on the container rather than 126 on the buttons. This is
        // event delegation, and it works because keyboard events BUBBLE up from
        // the focused button to its ancestors.
        on: { keydown: (event) => this.#handleKeyDown(event, ticks) },
      },
      buttons,
    );

    render(this.#root, [
      this.#rail,
      el(
        'p',
        { class: 'timeline__hint' },
        'Bar height shows how many events are recorded in that year. Use the arrow keys to ' +
          'move year by year, Page Up and Page Down to jump a decade.',
      ),
    ]);

    this.#renderedTicks = ticks;
  }

  /** Changes only the buttons that differ — two elements, not 126. */
  #applySelection(year) {
    const previous = this.#buttons.get(this.#selectedYear);
    if (previous) {
      previous.classList.remove('timeline__year--selected');
      previous.setAttribute('aria-checked', 'false');
      previous.setAttribute('tabindex', '-1');
    }

    const next = this.#buttons.get(year);
    if (next) {
      next.classList.add('timeline__year--selected');
      next.setAttribute('aria-checked', 'true');
      // ROVING TABINDEX: exactly one year is in the tab order and the arrow keys
      // move between them. Tabbing through 126 buttons is not navigation, it is
      // a punishment.
      next.setAttribute('tabindex', '0');
      next.scrollIntoView({ behavior: 'smooth', block: 'nearest', inline: 'center' });

      // Only move focus if focus is already inside the rail. Stealing focus from
      // elsewhere on the page whenever data loads would be hostile.
      if (this.#rail?.contains(document.activeElement) && document.activeElement !== next) {
        next.focus();
      }
    }

    this.#selectedYear = year;
  }

  #handleKeyDown(event, ticks) {
    if (this.#selectedYear === null || ticks.length === 0) return;

    const first = ticks[0].year;
    const last = ticks[ticks.length - 1].year;

    const step = { ArrowRight: 1, ArrowLeft: -1, ArrowUp: -1, ArrowDown: 1 }[event.key];
    const jump = { PageUp: -10, PageDown: 10 }[event.key];

    let next = null;
    if (step !== undefined) next = this.#selectedYear + step;
    else if (jump !== undefined) next = this.#selectedYear + jump;
    else if (event.key === 'Home') next = first;
    else if (event.key === 'End') next = last;
    if (next === null) return;

    // Stop the arrow keys from also scrolling the page.
    event.preventDefault();
    // Clamp instead of running off the end of the timeline.
    this.#onSelectYear(Math.min(last, Math.max(first, next)));
  }
}
