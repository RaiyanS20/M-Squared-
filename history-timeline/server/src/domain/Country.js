import { ValidationError } from './errors.js';

/**
 * A country users can browse.
 *
 * OOP lesson: ENTITY.
 *
 * Unlike a `Year` (a value object), a Country has an IDENTITY: two countries are
 * the same country when their `id` matches, even if someone later fixes a typo
 * in the name. Entity vs value object is the most useful distinction in domain
 * modelling, and getting it backwards causes real bugs — compare two entities
 * field-by-field and a renamed country looks like a *different* country.
 *
 * NOTE the different privacy technique from `Year`. Here the fields are public
 * and the whole object is frozen. `Object.freeze` is less strict than `#private`
 * (you can still read the fields, and freezing is shallow) but it is far less
 * boilerplate for an object with four fields, and it gives us what we need:
 * nobody can change a Country after it is built. Choosing the lighter tool when
 * it is sufficient is a real engineering skill; both approaches are shown here
 * on purpose.
 */
export class Country {
  constructor({ id, code, name, region }) {
    if (typeof name !== 'string' || name.trim() === '') {
      throw new ValidationError('Country name is required.', 'name');
    }
    if (typeof code !== 'string' || !/^[A-Za-z]{3}$/.test(code)) {
      throw new ValidationError(
        `Country code must be 3 letters (ISO 3166-1 alpha-3), received ${JSON.stringify(code)}.`,
        'code',
      );
    }

    this.id = id;
    this.code = code.toUpperCase();
    this.name = name.trim();
    this.region = typeof region === 'string' && region.trim() ? region.trim() : 'Unknown';

    // Freeze AFTER assigning. In a subclass you would freeze in the most
    // derived constructor instead, or the parent would lock the object before
    // the child could finish building it.
    Object.freeze(this);
  }

  /** Entities compare by identity. */
  equals(other) {
    return other instanceof Country && other.id === this.id;
  }

  /**
   * The wire format, written out explicitly.
   *
   * Letting JSON.stringify dump every field means the day someone adds an
   * internal note to this class, it silently ships to every browser. An explicit
   * toJSON is a deliberate boundary between your domain and the outside world.
   */
  toJSON() {
    return { id: this.id, code: this.code, name: this.name, region: this.region };
  }
}
