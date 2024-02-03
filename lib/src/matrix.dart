import 'dart:typed_data';

import 'package:collection/collection.dart';
import 'package:lapack/src/box.dart';

class Array<T> {
  final int _offset;
  final List<T> _array;

  Array(int length, {int offset = 0})
      : _array = switch (T) {
          double => Float64List(length),
          int => Int64List(length),
          bool => List.filled(length, false),
          _ => throw UnimplementedError(),
        } as List<T>,
        _offset = offset;

  Array.fromList(List<T> list, {int offset = 0})
      : _offset = offset,
        _array = switch (T) {
          double => Float64List.fromList(list as List<double>),
          int => Int64List.fromList(list as List<int>),
          bool => [...list],
          _ => throw UnimplementedError(),
        } as List<T>;

  Array.fromSlice(this._array, {int offset = 0}) : _offset = offset;

  Array<T> slice(int index) {
    return Array.fromSlice(switch (T) {
      double => (_array as Float64List).buffer.asFloat64List(
          (_array as Float64List).offsetInBytes +
              (index - _offset) * (_array as Float64List).elementSizeInBytes),
      int => (_array as Int64List).buffer.asInt64List(
          (_array as Int64List).offsetInBytes +
              (index - _offset) * (_array as Int64List).elementSizeInBytes),
      bool => _array is ListSlice
          ? _array.slice(index - _offset, _array.length)
          : ListSlice(_array, index - _offset, _array.length),
      _ => throw UnimplementedError(),
    } as List<T>);
  }

  Array<T> call(int index) {
    return slice(index);
  }

  Array<T> oneIndexed({int offset = 1}) {
    return Array.fromSlice(_array, offset: offset);
  }

  T operator [](int index) {
    return _array[index - _offset];
  }

  void operator []=(int index, T value) {
    _array[index - _offset] = value;
  }

  ArrayBox<T> box(index) {
    return ArrayBox(this, index);
  }
}

class Matrix2d<T> extends Array<T> {
  Matrix2d(int m, int n) : super(m * n);

  Matrix2d.fromList(List<T> list, {int offset = 0})
      : super.fromList(list, offset: offset);

  Matrix2d.fromArray(Array<T> array, {int offset = 0})
      : super.fromSlice(array._array, offset: offset);

  @override
  Matrix2d<T> slice(int index) {
    return Matrix2d.fromArray(super.slice(index));
  }

  @override
  Matrix2d<T> oneIndexed({int offset = 1}) {
    return Matrix2d.fromArray(this, offset: offset);
  }
}

class ArrayBox<T> implements Box<T> {
  final Array<T> _value;
  final int _index;
  ArrayBox(this._value, this._index);

  @override
  T get value => _value[_index];
  @override
  set value(T value) => _value[_index] = value;

  @override
  String toString() => value.toString();
}
