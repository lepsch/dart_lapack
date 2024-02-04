import 'package:complex/complex.dart' as impl;

class Complex {
  static const zero = Complex(0.0);
  static const one = Complex(1.0);

  final double real;
  final double imaginary;

  const Complex(this.real, [this.imaginary = 0]);

  factory Complex._fromImpl(impl.Complex impl) {
    return Complex(impl.real, impl.imaginary);
  }

  Complex operator +(Complex other) {
    return Complex._fromImpl(_impl + other._impl);
  }

  Complex operator -(Complex other) {
    return Complex._fromImpl(_impl - other._impl);
  }

  Complex operator *(Complex other) {
    return Complex._fromImpl(_impl * other._impl);
  }

  Complex operator /(Complex other) {
    return Complex._fromImpl(_impl / other._impl);
  }

  @override
  bool operator ==(Object? other) =>
      identical(this, other) ||
      other is Complex && real == other.real && imaginary == other.imaginary;

  @override
  int get hashCode => Object.hash(real, imaginary);

  impl.Complex get _impl => impl.Complex(real, imaginary);

  Complex conjugate() => Complex._fromImpl(_impl.conjugate());

  double toDouble() => real;
}

extension DoubleComplexExtension on double {
  Complex toComplex() => Complex(this);
}
