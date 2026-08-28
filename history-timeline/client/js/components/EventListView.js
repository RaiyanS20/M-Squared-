import { el, render } from '../dom.js';
import { loadingMessage, errorMessage, emptyMessage } from './statusMessages.js';

/**
 * Turns "04-27" into "27 April".
 *
 * A pure function — same input, same output, no DOM, no state. Pure functions
 * are the easiest thing in software to test, which is why this is exported
 * separately rather than buried inside the class.
 */
export function formatMonthDay(monthDay) {
  if (!monthDay) return null;
  const [month, day] = String(monthDay).split('-').map(Number);
  if (!month || !day || month < 1 || month > 12 || day < 1 || day > 31) return null;
  const monthName = new Date(Date.UTC(2000, month - 1, 1)).toLocaleString('en-GB', {
    month: 'long',
    timeZone: 'UTC',
  });
  return `${day} ${monthName}`;
}

/** A machine-readable value for <time>: "1994-04-27", or just "1994". */
function dateTimeAttr(year, monthDay) {
  return monthDay ? `${year}-${monthDay}` : String(year);
}

/**
 * MVP feature 4: the payoff.
 *
 * Note the source link on every card. That is not decoration — it is the whole
 * editorial promise of the project made visible, and the domain layer on the
 * server refuses to create an event that could not render one.
 */
export class EventListView {
  #root;

  constructor(root) {
    this.#root = root;
  }

  update(state) {
    const { events } = state;

    if (events.status === 'idle') {
      render(this.#root, emptyMessage('Choose a country to see its events.'));
      return;
    }
    if (events.status === 'loading') {
      render(this.#root, loadingMessage('Loading events…'));
      return;
    }
    if (events.status === 'error') {
      render(this.#root, errorMessage(events.error));
      return;
    }

    const { country, year, events: list } = events.data;

    if (list.length === 0) {
      render(this.#root, emptyMessage(`Nothing is recorded for ${country.name} in ${year}.`));
      return;
    }

    render(
      this.#root,
      el(
        'ol',
        { class: 'events', 'aria-label': `Events in ${country.name} in ${year}` },
        list.map((event) => {
          const date = formatMonthDay(event.monthDay);
          return el('li', { class: 'events__item' }, [
            el('article', { class: 'event' }, [
              el('header', { class: 'event__header' }, [
                el('span', { class: `event__category event__category--${event.category}` }, event.category),
                el(
                  'time',
                  { class: 'event__date', dateTime: dateTimeAttr(event.year, event.monthDay) },
                  date ? `${date} ${event.year}` : String(event.year),
                ),
              ]),
              el('h3', { class: 'event__title' }, event.title),
              el('p', { class: 'event__summary' }, event.summary),
              el(
                'a',
                {
                  class: 'event__source',
                  href: event.sourceUrl,
                  target: '_blank',
                  // `noopener` stops the opened page reaching back through
                  // window.opener; `noreferrer` withholds the referring URL.
                  // Always use both on target="_blank".
                  rel: 'noopener noreferrer',
                },
                [
                  'Read the source',
                  el(
                    'span',
                    { class: 'visually-hidden' },
                    ` for ${event.title} (opens in a new tab)`,
                  ),
                ],
              ),
            ]),
          ]);
        }),
      ),
    );
  }
}
