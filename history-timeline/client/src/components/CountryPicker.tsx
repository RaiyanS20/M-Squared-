import type { Country } from '../api/index.js';
import { useCountriesByRegion } from '../hooks/index.js';

interface CountryPickerProps {
  countries: Country[];
  selectedCode: string | null;
  year: number;
  onSelectCountry: (code: string) => void;
}

/**
 * MVP feature 3: list the countries, grouped by region.
 *
 * Only countries with something recorded in the chosen year reach this
 * component — the API already filtered them — so every option leads somewhere.
 * That decision was made in the SERVICE layer, not here; the UI just renders
 * what it is given.
 */
export function CountryPicker({
  countries,
  selectedCode,
  year,
  onSelectCountry,
}: CountryPickerProps) {
  const byRegion = useCountriesByRegion(countries);

  return (
    <nav className="countries" aria-label={`Countries with recorded events in ${year}`}>
      {byRegion.map(([region, list]) => (
        <section key={region} className="countries__region">
          <h3 className="countries__region-name">{region}</h3>
          <ul className="countries__list">
            {list.map((c) => {
              const isSelected = c.code === selectedCode;
              return (
                <li key={c.code}>
                  <button
                    type="button"
                    className={`countries__item${isSelected ? ' countries__item--selected' : ''}`}
                    aria-current={isSelected ? 'true' : undefined}
                    onClick={() => onSelectCountry(c.code)}
                  >
                    <span className="countries__name">{c.name}</span>
                    <span className="countries__code">{c.code}</span>
                  </button>
                </li>
              );
            })}
          </ul>
        </section>
      ))}
    </nav>
  );
}
