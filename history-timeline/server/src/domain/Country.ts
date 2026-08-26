import { ValidationError } from './errors.js';

/** The raw shape a Country arrives in from the database or a seed file. */
export interface CountryProps {
  id: number;
  /** ISO 3166-1 alpha-3 code, e.g. "FRA". Uppercase. */
  code: string;
  name: string;
  /** UN geoscheme region, e.g. "Europe". Used to group the country list. */
  region: string;
}

/**
 * A country users can browse.
 *
 * OOP lesson: ENTITY.
 * Unlike a `Year` (a value object), a Country has an IDENTITY: two countries are
 * the same country if their `id` matches, even if someone later fixes a typo in
 * the name. That is the difference between an entity and a value object, and it
 * is the single most useful distinction in domain modelling.
 */
export class Country {
  readonly id: number;
  readonly code: string;
  readonly name: string;
  readonly region: string;

  constructor(props: CountryProps) {
    if (!props.name?.trim()) {
      throw new ValidationError('Country name is required.', 'name');
    }
    if (!/^[A-Za-z]{3}$/.test(props.code ?? '')) {
      throw new ValidationError(
        `Country code must be 3 letters (ISO 3166-1 alpha-3), received "${props.code}".`,
        'code',
      );
    }
    this.id = props.id;
    this.code = props.code.toUpperCase();
    this.name = props.name.trim();
    this.region = props.region?.trim() || 'Unknown';
    Object.freeze(this);
  }

  /** Entities compare by identity. */
  equals(other: Country): boolean {
    return other instanceof Country && other.id === this.id;
  }

  /**
   * The wire format. Keeping this explicit (rather than letting `JSON.stringify`
   * dump every field) means a private field added later is not accidentally
   * leaked to the browser. This is the boundary between your domain and your API.
   */
  toJSON() {
    return { id: this.id, code: this.code, name: this.name, region: this.region };
  }
}
