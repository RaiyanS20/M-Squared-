import { describe, expect, it, vi } from 'vitest';
import { render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { Timeline } from '../src/components/Timeline.js';
import { CountryPicker } from '../src/components/CountryPicker.js';
import { EventList, formatMonthDay } from '../src/components/EventList.js';
import { EVENTS_FRA_1969, FRANCE, JAPAN, TIMELINE, anEvent } from './stubApi.js';

/**
 * Component tests written the way a USER experiences the component: find things
 * by their visible text or accessible role, click them, assert on what appears.
 *
 * Nothing here asserts on a CSS class or a component's internal state. That is
 * deliberate — those are implementation details, and tests coupled to them break
 * on every refactor while catching none of the bugs that matter.
 */

describe('Timeline', () => {
  it('renders one radio per year with an accessible label', () => {
    render(<Timeline ticks={TIMELINE.ticks} selectedYear={1969} onSelectYear={vi.fn()} />);
    expect(screen.getAllByRole('radio')).toHaveLength(4);
    expect(screen.getByRole('radio', { name: /1969, 2 events/ })).toBeInTheDocument();
    expect(screen.getByRole('radio', { name: /1968, 1 event$/ })).toBeInTheDocument();
  });

  it('marks the selected year as checked, and no other', () => {
    render(<Timeline ticks={TIMELINE.ticks} selectedYear={1969} onSelectYear={vi.fn()} />);
    const checked = screen.getAllByRole('radio').filter((r) => r.getAttribute('aria-checked') === 'true');
    expect(checked).toHaveLength(1);
    expect(checked[0]).toHaveAccessibleName(/1969/);
  });

  it('reports the year that was clicked', async () => {
    const onSelectYear = vi.fn();
    render(<Timeline ticks={TIMELINE.ticks} selectedYear={1969} onSelectYear={onSelectYear} />);

    await userEvent.click(screen.getByRole('radio', { name: /1967/ }));
    expect(onSelectYear).toHaveBeenCalledWith(1967);
  });

  describe('keyboard navigation', () => {
    it('moves one year with the arrow keys', async () => {
      const onSelectYear = vi.fn();
      render(<Timeline ticks={TIMELINE.ticks} selectedYear={1969} onSelectYear={onSelectYear} />);

      screen.getByRole('radio', { name: /1969/ }).focus();
      await userEvent.keyboard('{ArrowRight}');
      expect(onSelectYear).toHaveBeenCalledWith(1970);

      await userEvent.keyboard('{ArrowLeft}');
      expect(onSelectYear).toHaveBeenLastCalledWith(1968);
    });

    it('clamps at the ends of the range instead of running off', async () => {
      const onSelectYear = vi.fn();
      render(<Timeline ticks={TIMELINE.ticks} selectedYear={1970} onSelectYear={onSelectYear} />);

      screen.getByRole('radio', { name: /1970/ }).focus();
      await userEvent.keyboard('{ArrowRight}');
      expect(onSelectYear).toHaveBeenCalledWith(1970);
    });

    it('jumps to the first and last year with Home and End', async () => {
      const onSelectYear = vi.fn();
      render(<Timeline ticks={TIMELINE.ticks} selectedYear={1969} onSelectYear={onSelectYear} />);

      screen.getByRole('radio', { name: /1969/ }).focus();
      await userEvent.keyboard('{Home}');
      expect(onSelectYear).toHaveBeenLastCalledWith(1967);
      await userEvent.keyboard('{End}');
      expect(onSelectYear).toHaveBeenLastCalledWith(1970);
    });

    it('keeps exactly one year in the tab order', () => {
      render(<Timeline ticks={TIMELINE.ticks} selectedYear={1969} onSelectYear={vi.fn()} />);
      const tabbable = screen.getAllByRole('radio').filter((r) => r.getAttribute('tabindex') === '0');
      expect(tabbable).toHaveLength(1);
    });
  });
});

describe('CountryPicker', () => {
  it('groups countries under their region', () => {
    render(
      <CountryPicker countries={[FRANCE, JAPAN]} selectedCode={null} year={1969} onSelectCountry={vi.fn()} />,
    );
    expect(screen.getByRole('heading', { name: 'Asia' })).toBeInTheDocument();
    expect(screen.getByRole('heading', { name: 'Europe' })).toBeInTheDocument();
    expect(screen.getByRole('button', { name: /France/ })).toBeInTheDocument();
  });

  it('reports the chosen country code', async () => {
    const onSelectCountry = vi.fn();
    render(
      <CountryPicker countries={[FRANCE, JAPAN]} selectedCode={null} year={1969} onSelectCountry={onSelectCountry} />,
    );
    await userEvent.click(screen.getByRole('button', { name: /Japan/ }));
    expect(onSelectCountry).toHaveBeenCalledWith('JPN');
  });

  it('marks the selected country for assistive technology', () => {
    render(
      <CountryPicker countries={[FRANCE, JAPAN]} selectedCode="FRA" year={1969} onSelectCountry={vi.fn()} />,
    );
    expect(screen.getByRole('button', { name: /France/ })).toHaveAttribute('aria-current', 'true');
  });
});

describe('formatMonthDay', () => {
  it.each([
    ['03-02', '2 March'],
    ['12-25', '25 December'],
    ['01-01', '1 January'],
  ])('formats %s as %s', (input, expected) => {
    expect(formatMonthDay(input)).toBe(expected);
  });

  it.each([null, '', '13-01', 'nonsense'])('returns null for %o', (input) => {
    expect(formatMonthDay(input as string | null)).toBeNull();
  });
});

describe('EventList', () => {
  it('renders each event with its date, summary and category', () => {
    render(<EventList events={EVENTS_FRA_1969.events} countryName="France" year={1969} />);
    expect(screen.getByText('Concorde first flight')).toBeInTheDocument();
    expect(screen.getByText('2 March 1969')).toBeInTheDocument();
    expect(screen.getByText('science')).toBeInTheDocument();
    expect(screen.getAllByRole('listitem')).toHaveLength(2);
  });

  it('shows just the year when the exact date is unknown', () => {
    render(<EventList events={[anEvent({ monthDay: null })]} countryName="France" year={1969} />);
    expect(screen.getByText('1969')).toBeInTheDocument();
  });

  // The editorial promise, verified in the UI.
  it('gives every event a source link that opens safely in a new tab', () => {
    render(<EventList events={EVENTS_FRA_1969.events} countryName="France" year={1969} />);
    const links = screen.getAllByRole('link');
    expect(links).toHaveLength(EVENTS_FRA_1969.events.length);
    for (const link of links) {
      expect(link).toHaveAttribute('href', expect.stringMatching(/^https?:\/\//));
      expect(link).toHaveAttribute('rel', 'noopener noreferrer');
    }
  });

  it('exposes a machine-readable date', () => {
    const { container } = render(
      <EventList events={[anEvent({ monthDay: '04-28', year: 1969 })]} countryName="France" year={1969} />,
    );
    expect(container.querySelector('time')).toHaveAttribute('datetime', '1969-04-28');
  });
});
