import { describe, expect, it } from 'vitest';
import { render, screen, waitFor, within } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { App } from '../src/App.js';
import { COUNTRIES_1969, EVENTS_FRA_1969, TIMELINE, stubFetch } from './stubApi.js';

/**
 * END-TO-END-ISH: the whole app, driven the way a learner uses it, with only the
 * network stubbed. This one file proves all four MVP features work together —
 * which is the assertion that actually matters to the product.
 */
describe('App', () => {
  const happyPath = () =>
    stubFetch({
      '/api/timeline': TIMELINE,
      '/api/years/1969/countries/FRA/events': EVENTS_FRA_1969,
      '/api/years/1969/countries': COUNTRIES_1969,
    });

  it('walks a learner from year to country to events', async () => {
    happyPath();
    render(<App />);

    // Feature 1: the timeline appears, and the app pre-selects the most recent
    // year that actually has events (1969, not the empty 1970).
    const selected = await screen.findByRole('radio', { checked: true });
    expect(selected).toHaveAccessibleName(/1969/);

    // Feature 3: countries for that year.
    expect(await screen.findByRole('button', { name: /France/ })).toBeInTheDocument();
    expect(screen.getByRole('button', { name: /Japan/ })).toBeInTheDocument();

    // Feature 4: pick one and read what happened.
    await userEvent.click(screen.getByRole('button', { name: /France/ }));
    expect(await screen.findByText('Concorde first flight')).toBeInTheDocument();
    expect(screen.getByText('De Gaulle resigns')).toBeInTheDocument();
  });

  it('prompts for a country before any is chosen', async () => {
    happyPath();
    render(<App />);
    expect(await screen.findByText(/choose a country to see its events/i)).toBeInTheDocument();
  });

  it('clears the chosen country when the year changes', async () => {
    happyPath();
    render(<App />);

    await userEvent.click(await screen.findByRole('button', { name: /France/ }));
    expect(await screen.findByText('Concorde first flight')).toBeInTheDocument();

    // Feature 2: move to a different year. The previous country's events must
    // not linger — that is the stale-state bug this test exists to prevent.
    await userEvent.click(screen.getByRole('radio', { name: /1967/ }));

    await waitFor(() => {
      expect(screen.queryByText('Concorde first flight')).not.toBeInTheDocument();
    });
    expect(screen.getByText(/choose a country to see its events/i)).toBeInTheDocument();
  });

  it('explains a year with nothing recorded instead of showing a blank panel', async () => {
    stubFetch({
      '/api/timeline': TIMELINE,
      '/api/years/1967/countries': { year: 1967, countries: [] },
      '/api/years/1969/countries': COUNTRIES_1969,
    });
    render(<App />);

    await userEvent.click(await screen.findByRole('radio', { name: /1967/ }));
    expect(await screen.findByText(/nothing is recorded for 1967/i)).toBeInTheDocument();
  });

  it('surfaces a server error instead of failing silently', async () => {
    stubFetch(
      { '/api/timeline': { error: { code: 'INTERNAL_ERROR', message: 'The archive is offline.' } } },
      { status: { '/api/timeline': 500 } },
    );
    render(<App />);

    const alert = await screen.findByRole('alert');
    expect(within(alert).getByText('The archive is offline.')).toBeInTheDocument();
  });

  it('has one labelled step per MVP feature', async () => {
    happyPath();
    render(<App />);
    expect(screen.getByRole('heading', { name: /1\. choose a year/i })).toBeInTheDocument();
    expect(screen.getByRole('heading', { name: /2\. choose a country/i })).toBeInTheDocument();
    expect(screen.getByRole('heading', { name: /3\. what happened/i })).toBeInTheDocument();
  });
});
