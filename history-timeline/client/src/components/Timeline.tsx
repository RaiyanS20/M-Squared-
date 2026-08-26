import { useEffect, useRef } from 'react';
import type { TimelineTick } from '../api/index.js';

interface TimelineProps {
  ticks: TimelineTick[];
  selectedYear: number | null;
  onSelectYear: (year: number) => void;
}

/**
 * MVP features 1 and 2: the scale, and choosing a year on it.
 *
 * Three things here are worth copying into any project:
 *
 *  1. It is a CONTROLLED component. It holds no state of its own — the selected
 *     year is passed in, and a click calls back up. State lives in one place
 *     (App), so the timeline, the country list and the URL can never disagree.
 *
 *  2. It is a real ARIA radiogroup. A row of divs with onClick is invisible to
 *     assistive technology and unusable from a keyboard. Native buttons plus
 *     roving tabindex give arrow-key navigation for free.
 *
 *  3. The density bar is derived, never stored. `eventCount` comes from the API;
 *     the bar height is computed at render time from the busiest year in view.
 */
export function Timeline({ ticks, selectedYear, onSelectYear }: TimelineProps) {
  const railRef = useRef<HTMLDivElement>(null);
  const busiest = Math.max(1, ...ticks.map((t) => t.eventCount));

  // Keep the chosen year in view — including on first load, where the default
  // selection may be hundreds of pixels off-screen.
  useEffect(() => {
    if (selectedYear === null) return;
    railRef.current
      ?.querySelector(`[data-year="${selectedYear}"]`)
      ?.scrollIntoView({ behavior: 'smooth', block: 'nearest', inline: 'center' });
  }, [selectedYear]);

  const handleKeyDown = (event: React.KeyboardEvent<HTMLDivElement>) => {
    if (selectedYear === null) return;
    const step = { ArrowRight: 1, ArrowLeft: -1, ArrowUp: -1, ArrowDown: 1 }[event.key];
    const jump = { PageUp: -10, PageDown: 10 }[event.key];
    const first = ticks[0]?.year;
    const last = ticks.at(-1)?.year;
    if (first === undefined || last === undefined) return;

    let next: number | null = null;
    if (step !== undefined) next = selectedYear + step;
    else if (jump !== undefined) next = selectedYear + jump;
    else if (event.key === 'Home') next = first;
    else if (event.key === 'End') next = last;
    if (next === null) return;

    event.preventDefault();
    onSelectYear(Math.min(last, Math.max(first, next)));
  };

  return (
    <div className="timeline">
      <div
        className="timeline__rail"
        ref={railRef}
        role="radiogroup"
        aria-label="Choose a year"
        onKeyDown={handleKeyDown}
      >
        {ticks.map((tick) => {
          const isSelected = tick.year === selectedYear;
          const isDecade = tick.year % 10 === 0;
          return (
            <button
              key={tick.year}
              type="button"
              role="radio"
              aria-checked={isSelected}
              // Roving tabindex: exactly one year is in the tab order, and the
              // arrow keys move between them. Tabbing through 127 buttons is not
              // navigation, it is a punishment.
              tabIndex={isSelected ? 0 : -1}
              data-year={tick.year}
              className={[
                'timeline__year',
                isSelected && 'timeline__year--selected',
                isDecade && 'timeline__year--decade',
                tick.eventCount === 0 && 'timeline__year--empty',
              ]
                .filter(Boolean)
                .join(' ')}
              onClick={() => onSelectYear(tick.year)}
              aria-label={`${tick.year}, ${tick.eventCount} ${
                tick.eventCount === 1 ? 'event' : 'events'
              }`}
              title={`${tick.year} — ${tick.eventCount} recorded ${
                tick.eventCount === 1 ? 'event' : 'events'
              }`}
            >
              <span
                className="timeline__bar"
                style={{ height: `${(tick.eventCount / busiest) * 100}%` }}
                aria-hidden="true"
              />
              <span className="timeline__label" aria-hidden="true">
                {isDecade ? tick.year : ''}
              </span>
            </button>
          );
        })}
      </div>
      <p className="timeline__hint">
        Bar height shows how many events are recorded in that year. Use the arrow keys to move year
        by year, Page Up and Page Down to jump a decade.
      </p>
    </div>
  );
}
