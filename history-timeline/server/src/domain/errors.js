/**
 * Domain errors.
 *
 * OOP lesson: INHERITANCE + POLYMORPHISM.
 *
 * Every error our business rules can raise extends `DomainError`. That lets one
 * `catch` in the HTTP layer ask "is this a domain problem, or a bug?" without
 * knowing every concrete subclass. Adding a new error type later needs zero
 * changes to the error handler.
 */
export class DomainError extends Error {
  constructor(message) {
    super(message);

    // JavaScript has no `abstract` keyword, so we enforce it at runtime:
    // `new DomainError(...)` is a mistake — every error should be specific.
    if (new.target === DomainError) {
      throw new TypeError('DomainError is abstract; throw a subclass instead.');
    }

    // `new.target` is the constructor that was actually called, so an instance
    // of ValidationError reports `name === 'ValidationError'`.
    this.name = new.target.name;
    Error.captureStackTrace?.(this, new.target);
  }
}

/** The caller sent something that cannot become a valid domain object. */
export class ValidationError extends DomainError {
  /** Machine-readable code. The client switches on this, never on the message. */
  code = 'VALIDATION_ERROR';
  /** The status the API layer should map this to. */
  httpStatus = 400;

  constructor(message, field) {
    super(message);
    this.field = field;
  }
}

/** The caller asked for something that does not exist. */
export class NotFoundError extends DomainError {
  code = 'NOT_FOUND';
  httpStatus = 404;

  constructor(resource, identifier) {
    super(`${resource} '${identifier}' was not found.`);
  }
}
