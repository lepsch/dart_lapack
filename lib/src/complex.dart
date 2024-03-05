import 'dart:collection';
import 'dart:typed_data';

import 'package:complex/complex.dart' as impl;

abstract class Complex {
  static const zero = Complex(0.0);
  static const one = Complex(1.0);
  static const nan = Complex(double.nan, double.nan);

  const factory Complex(double real, [double imaginary]) = _Complex;

  factory Complex.fromFloat64x2(Float64x2 float) {
    return Complex(float.x, float.y);
  }

  double get real;
  double get imaginary;
  bool get isNaN;

  Complex operator +(Complex other);
  Complex operator -(Complex other);
  Complex operator -();
  Complex operator *(Complex other);
  Complex operator /(Complex other);
  bool operator <(Complex other);
  bool operator <=(Complex other);
  bool operator >(Complex other);
  bool operator >=(Complex other);
  @override
  bool operator ==(Object? other);
  @override
  int get hashCode;

  Complex conjugate();

  Complex sqrt();

  Complex pow(num x);

  double toDouble() => real;

  int toInt() => real.toInt();

  double abs();

  Complex exp();

  Complex sin();

  Complex cos();

  @override
  String toString();
}

class _Complex implements Complex {
  @override
  final double real;
  @override
  final double imaginary;
  @override
  bool get isNaN => real.isNaN && imaginary.isNaN;

  const _Complex(this.real, [this.imaginary = 0]);

  factory _Complex._fromImpl(impl.Complex impl) {
    return _Complex(impl.real, impl.imaginary);
  }

  @override
  Complex operator +(Complex other) {
    return _Complex._fromImpl(
        _impl + impl.Complex(other.real, other.imaginary));
  }

  @override
  Complex operator -(Complex other) {
    return _Complex._fromImpl(
        _impl - impl.Complex(other.real, other.imaginary));
  }

  @override
  Complex operator -() {
    return _Complex._fromImpl(-_impl);
  }

  @override
  Complex operator *(Complex other) {
    return _Complex._fromImpl(
        _impl * impl.Complex(other.real, other.imaginary));
  }

  @override
  Complex operator /(Complex other) {
    return _Complex._fromImpl(
        _impl / impl.Complex(other.real, other.imaginary));
  }

  @override
  bool operator <(Complex other) {
    return abs() < other.abs();
  }

  @override
  bool operator <=(Complex other) {
    return abs() <= other.abs();
  }

  @override
  bool operator >(Complex other) {
    return abs() > other.abs();
  }

  @override
  bool operator >=(Complex other) {
    return abs() >= other.abs();
  }

  @override
  bool operator ==(Object? other) =>
      identical(this, other) ||
      other is Complex && real == other.real && imaginary == other.imaginary ||
      other is ComplexTuple && real == other.$1 && imaginary == other.$2 ||
      other is num && real == other && imaginary == 0.0;

  @override
  int get hashCode => Object.hash(real, imaginary);

  impl.Complex get _impl => impl.Complex(real, imaginary);

  @override
  Complex conjugate() => _Complex._fromImpl(_impl.conjugate());

  @override
  Complex sqrt() => _Complex._fromImpl(_impl.sqrt());

  @override
  Complex pow(num x) {
    if (x is int) {
      if (this == Complex.zero) {
        return x == 0 ? Complex.one : Complex.zero;
      }
      if (x == 1) return this;

      // final result = _Complex._fromImpl(_impl.pow(x));
      // return Complex(
      //   result.real.toInt().toDouble(),
      //   result.imaginary.toInt().toDouble(),
      // );
    } else {
      if (this == Complex.zero) {
        if (x != 0.0) return Complex(0, -0.0);
      }
    }
    return _Complex._fromImpl(_impl.pow(x));
  }

  @override
  double toDouble() => real;

  @override
  int toInt() => real.toInt();

  @override
  double abs() => _impl.abs();

  @override
  Complex exp() => _Complex._fromImpl(_impl.exp());

  @override
  Complex sin() => _Complex._fromImpl(_impl.sin());

  @override
  Complex cos() => _Complex._fromImpl(_impl.cos());

  @override
  String toString() {
    return '(${real.toStringAsExponential(4)}, ${imaginary.toStringAsExponential(4)})';
  }
}

class Complex64List with ListMixin<Complex> implements List<Complex> {
  final Float64x2List _list;

  Complex64List(int length) : _list = Float64x2List(length);
  Complex64List._(this._list);

  factory Complex64List.fromList(List<Complex> list) {
    if (list is Complex64List) {
      return Complex64List._(Float64x2List.fromList(list._list));
    }
    return Complex64List._(
        Float64x2List.fromList(list.map((c) => c.toFloat64x2()).toList()));
  }

  Complex64List.fromData(this._list);

  factory Complex64List.fromFloat64List(List<Float64x2> list) {
    return Complex64List._(Float64x2List.fromList(list));
  }

  factory Complex64List.fromTuples(List<(double x, double y)> list) {
    return Complex64List._(Float64x2List.fromList(
        list.map((c) => Float64x2(c.$1, c.$2)).toList()));
  }

  Complex64List slice(int start, [int? length]) {
    return Complex64List._(_list.buffer.asFloat64x2List(
        _list.offsetInBytes + start * _list.elementSizeInBytes, length));
  }

  Float64x2List toData() => _list;

  @override
  int get length => _list.length;

  @override
  set length(int newLength) {
    _list.length = newLength;
  }

  @override
  Complex operator [](int index) {
    return _list[index].toComplex();
  }

  @override
  void operator []=(int index, Complex value) {
    _list[index] = value.toFloat64x2();
  }
}

extension DoubleComplexExtension on double {
  Complex toComplex() => Complex(this);
}

extension NumComplexExtension on num {
  Complex toComplex() => Complex(toDouble());
}

extension Float64x2ComplexExtension on Float64x2 {
  Complex toComplex() => Complex.fromFloat64x2(this);
}

extension ComplexFloat64x2Extension on Complex {
  Float64x2 toFloat64x2() => Float64x2(real, imaginary);
}

typedef ComplexTuple = (num, num);

extension ComplexTupleExtension on ComplexTuple {
  Complex toComplex() => Complex($1.toDouble(), $2.toDouble());
}

extension ComplexListTupleExtension on List<ComplexTuple> {
  List<Complex> toComplexList() => map((c) => c.toComplex()).toList();
}
