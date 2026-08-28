import { Store, Async } from './Store.js';
import { ApiError } from './api/ApiClient.js';
import { TimelineView } from './components/TimelineView.js';
import { CountryListView } from './components/CountryListView.js';
import { EventListView } from './components/EventListView.js';

/**
 * The app shell: it owns the state, loads the data, and tells the views to
 * update. It touches no DOM directly — the views do that.
 *
 * ---------------------------------------------------------------------------
 * STATE DESIGN — the most important decision in any front-end application.
 * ---------------------------------------------------------------------------
 *
 * Only TWO things are really chosen by the user:
 *
 *     selectedYear, selectedCountryCode
 *
 * The three data slices (timeline, countries, events) are DERIVED from those:
 * when a selection changes, the matching request is re-issued. There is no
 * second copy of anything, so no two panels can disagree.
 *
 * The rule to carry with you: if you can compute it, do not store it.
 */
export class App {
  #api;
  #store;
  #views = [];
  /** One AbortController per data slice, so a new request cancels the old one. */
  #inFlight = { countries: null, events: null, timeline: null };

  constructor({ api, roots }) {
    this.#api = api;

    this.#store = new Store({
      selectedYear: null,
      selectedCountryCode: null,
      timeline: Async.idle(),
      countries: Async.idle(),
      events: Async.idle(),
    });

    this.#views = [
      new TimelineView(roots.timeline, { onSelectYear: (year) => this.selectYear(year) }),
      new CountryListView(roots.countries, {
        contextRoot: roots.countriesContext,
        onSelectCountry: (code) => this.selectCountry(code),
      }),
      new EventListView(roots.events),
    ];

    // Every state change re-runs every view. Each view then decides how little
    // work it needs to do — see the comment at the top of TimelineView.
    this.#store.subscribe((state) => this.#render(state));
  }

  /** For tests and debugging. */
  get state() {
    return this.#store.getState();
  }

  /** Loads the timeline and paints the first frame. */
  async start() {
    this.#render(this.#store.getState());
    await this.#loadTimeline();
  }

  #render(state) {
    for (const view of this.#views) view.update(state);
  }

  // -- Selection -------------------------------------------------------------

  selectYear(year) {
    if (year === this.#store.getState().selectedYear) return;

    // Changing the year MUST clear the country. "France in 1969" is a valid
    // pair, but after switching to 1848 France may have nothing recorded — and
    // leaving the old country selected would show one year's events under
    // another year's heading.
    //
    // Resetting dependent state when its parent changes is the cure for a whole
    // family of stale-UI bugs, and setting `events` back to idle here is the
    // other half of it: without that, the PREVIOUS country's events would stay
    // on screen, looking entirely plausible and being completely wrong.
    this.#cancel('events');
    this.#store.setState({
      selectedYear: year,
      selectedCountryCode: null,
      events: Async.idle(),
    });

    void this.#loadCountries(year);
  }

  selectCountry(code) {
    if (code === this.#store.getState().selectedCountryCode) return;
    this.#store.setState({ selectedCountryCode: code });
    void this.#loadEvents(this.#store.getState().selectedYear, code);
  }

  // -- Loading ---------------------------------------------------------------

  /**
   * Starts a request for one slice, cancelling any request already in flight for
   * that same slice.
   *
   * THE RACE CONDITION THIS PREVENTS: click 1914, then quickly click 1969. Two
   * requests are in flight. If 1914's response happens to arrive SECOND, it
   * overwrites 1969's, and the panel now shows 1914's countries under a 1969
   * heading. This is invisible on a fast local network and very visible on a
   * phone. Aborting the previous request makes the bug impossible rather than
   * merely unlikely.
   */
  async #load(slice, request) {
    this.#cancel(slice);
    const controller = new AbortController();
    this.#inFlight[slice] = controller;

    this.#store.setState({ [slice]: Async.loading() });

    try {
      const data = await request(controller.signal);
      if (controller.signal.aborted) return;
      this.#store.setState({ [slice]: Async.success(data) });
    } catch (error) {
      // We cancelled it on purpose; that is not a failure to report.
      if (controller.signal.aborted || error?.name === 'AbortError') return;
      this.#store.setState({
        [slice]: Async.error(
          error instanceof ApiError
            ? error
            : new ApiError(0, 'UNKNOWN_ERROR', 'Something unexpected went wrong.'),
        ),
      });
    } finally {
      if (this.#inFlight[slice] === controller) this.#inFlight[slice] = null;
    }
  }

  #cancel(slice) {
    this.#inFlight[slice]?.abort();
    this.#inFlight[slice] = null;
  }

  async #loadTimeline() {
    await this.#load('timeline', (signal) => this.#api.getTimeline(undefined, undefined, signal));

    // Start on the most recent year that actually has events, so a first-time
    // visitor sees content rather than an empty panel.
    const { timeline, selectedYear } = this.#store.getState();
    if (timeline.status === 'success' && selectedYear === null) {
      const populated = [...timeline.data.ticks].reverse().find((t) => t.eventCount > 0);
      this.selectYear(populated?.year ?? timeline.data.endYear);
    }
  }

  async #loadCountries(year) {
    await this.#load('countries', (signal) => this.#api.getCountriesForYear(year, signal));
  }

  async #loadEvents(year, code) {
    await this.#load('events', (signal) => this.#api.getEvents(year, code, signal));
  }
}
