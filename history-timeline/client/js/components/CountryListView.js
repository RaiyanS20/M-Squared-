import { el, render } from '../dom.js';
import { loadingMessage, errorMessage, emptyMessage } from './statusMessages.js';

/** Groups countries by region, sorting both the regions and the names within. */
export function groupByRegion(countries) {
  const byRegion = new Map();
  for (const country of countries) {
    const bucket = byRegion.get(country.region) ?? [];
    bucket.push(country);
    byRegion.set(country.region, bucket);
  }
  return [...byRegion.entries()]
    .map(([region, list]) => [region, [...list].sort((a, b) => a.name.localeCompare(b.name))])
    .sort(([a], [b]) => a.localeCompare(b));
}

/**
 * MVP feature 3: list the countries, grouped by region.
 *
 * Only countries with something recorded in the chosen year reach this view —
 * the API already filtered them, so every option leads somewhere. That decision
 * was made in the SERVICE layer on the server, not here. The view renders what
 * it is given.
 *
 * This one rebuilds its whole list on change, unlike TimelineView. That is a
 * deliberate, and correct, difference: the list is a dozen elements and its
 * contents change completely whenever the year changes, so the bookkeeping that
 * pays for itself in the timeline would be pure overhead here. Match the
 * technique to the problem.
 */
export class CountryListView {
  #root;
  #contextRoot;
  #onSelectCountry;

  constructor(root, { contextRoot, onSelectCountry }) {
    this.#root = root;
    this.#contextRoot = contextRoot;
    this.#onSelectCountry = onSelectCountry;
  }

  update(state) {
    const { countries, selectedYear, selectedCountryCode } = state;

    // The " in 1994" after the panel heading.
    this.#contextRoot.textContent = selectedYear === null ? '' : ` in ${selectedYear}`;

    if (countries.status === 'idle') {
      render(this.#root, emptyMessage('Select a year to begin.'));
      return;
    }
    if (countries.status === 'loading') {
      render(this.#root, loadingMessage('Finding countries…'));
      return;
    }
    if (countries.status === 'error') {
      render(this.#root, errorMessage(countries.error));
      return;
    }

    const list = countries.data.countries;
    if (list.length === 0) {
      render(
        this.#root,
        emptyMessage(
          `Nothing is recorded for ${countries.data.year} yet. Try a year with a taller bar on the timeline.`,
        ),
      );
      return;
    }

    render(
      this.#root,
      el(
        'nav',
        { 'aria-label': `Countries with recorded events in ${countries.data.year}` },
        groupByRegion(list).map(([region, inRegion]) =>
          el('section', { class: 'countries__region' }, [
            el('h3', { class: 'countries__region-name' }, region),
            el(
              'ul',
              { class: 'countries__list' },
              inRegion.map((country) => {
                const isSelected = country.code === selectedCountryCode;
                return el('li', {}, [
                  el(
                    'button',
                    {
                      type: 'button',
                      class: `countries__item${isSelected ? ' countries__item--selected' : ''}`,
                      // aria-current tells a screen reader which one is active.
                      'aria-current': isSelected ? 'true' : null,
                      on: { click: () => this.#onSelectCountry(country.code) },
                    },
                    [
                      el('span', { class: 'countries__name' }, country.name),
                      el('span', { class: 'countries__code' }, country.code),
                    ],
                  ),
                ]);
              }),
            ),
          ]),
        ),
      ),
    );
  }
}
