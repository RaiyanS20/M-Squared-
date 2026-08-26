/**
 * Domain errors.
 *
 * OOP lesson: INHERITANCE + POLYMORPHISM.
 * Every error our business rules can raise extends `DomainError`. That lets one
 * `catch` block in the HTTP layer ask "is this a domain problem or a bug?"
 * without knowing every concrete subclass. Adding a new error type later
 * requires zero changes to the error handler.
 */
export abstract class DomainError extends Error {
  /** HTTP-agnostic machine code, e.g. "VALIDATION_ERROR". */
  abstract readonly code: string;

  /** The status the API layer should map this to. Kept here for convenience. */
  abstract readonly httpStatus: number;

  protected constructor(message: string) {
    super(message);
    // Restores the prototype chain when compiling to ES5-ish targets and makes
    // `instanceof` reliable. Also gives errors their real class name.
    this.name = new.target.name;
    Error.captureStackTrace?.(this, new.target);
  }
}

/** The caller sent something that cannot be turned into a valid domain object. */
export class ValidationError extends DomainError {
  readonly code = 'VALIDATION_ERROR';
  readonly httpStatus = 400;

  constructor(message: string, readonly field?: string) {
    super(message);
  }
}

/** The caller asked for something that does not exist. */
export class NotFoundError extends DomainError {
  readonly code = 'NOT_FOUND';
  readonly httpStatus = 404;

  constructor(resource: string, identifier: string | number) {
    super(`${resource} '${identifier}' was not found.`);
  }
}
