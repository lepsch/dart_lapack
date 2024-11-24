// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'dart:math';
import 'dart:typed_data';

import 'package:collection/collection.dart';
import 'package:dart_lapack/src/box.dart';
import 'package:dart_lapack/src/complex.dart';

const oneIndexedArrayOffset = -1;
const zeroIndexedArrayOffset = 0;
const oneIndexedMatrixOffset = (x: -1, y: -1);
const zeroIndexedMatrixOffset = (x: 0, y: 0);
const oneIndexedMatrix3dOffset = (x: -1, y: -1, z: -1);
const zeroIndexedMatrix3dOffset = (x: 0, y: 0, z: 0);

abstract interface class Array<T> implements Box<T>, List<T> {
  factory Array(
    int length, {
    int offset = oneIndexedArrayOffset,
  }) {
    return _Array<T>(length, offset: offset);
  }

  factory Array.fromList(
    List<T> list, {
    int offset = oneIndexedArrayOffset,
  }) {
    return _Array.fromList(list, offset: offset);
  }

  factory Array.fromData(List<T> list, {int offset}) = Array.fromList;

  factory Array.fromSlice(
    List<T> elements, {
    int offset = oneIndexedArrayOffset,
  }) {
    return _Array.fromData(elements, offset: offset);
  }

  void assign(Array<T> array);

  Array<T> slice(int index, {int? offset});

  Array<T> call(int index, {int? offset});

  Box<T> box(int index);

  List<T> toData();

  Matrix<T> asMatrix([int ld]);

  Array<T> having({int? length, int? offset});

  Array<T> copy();

  @override
  Array<R> cast<R>();

  int get offset;
}

extension ArrayExtension<T> on Array<T> {
  T maxval(int start, int end) {
    switch (T) {
      case const (double):
      case const (int):
        var value = this[start] as num;
        for (var i = start + 1; i <= end; i++) {
          if (value < (this[i] as num)) value = this[i] as num;
        }
        return value as T;

      case const (Complex):
        var value = this[start] as Complex;
        for (var i = start + 1; i <= end; i++) {
          if (value < (this[i] as Complex)) value = this[i] as Complex;
        }
        return value as T;

      case const (bool):
        for (var i = start; i <= end; i++) {
          if (this[start] as bool) return true as T;
        }
        return false as T;
      default:
        throw UnimplementedError();
    }
  }

  int maxloc(int start, int end, {int dim = 1}) {
    switch (T) {
      case const (double):
      case const (int):
        int loc = start;
        var value = this[start] as num;
        for (var i = start + 1; i <= end; i++) {
          if (value < (this[i] as num)) {
            value = this[i] as num;
            loc = i;
          }
        }
        return loc - start + 1;

      case const (Complex):
        int loc = start;
        Complex value = this[start] as Complex;
        for (var i = start + 1; i <= end; i++) {
          if (value < (this[i] as Complex)) {
            value = this[i] as Complex;
            loc = i;
          }
        }
        return loc - start + 1;

      case const (bool):
        for (var i = start; i <= end; i++) {
          if (this[start] as bool) return i - start + 1;
        }
        return 0;
      default:
        throw UnimplementedError();
    }
  }
}

class Matrix<T> implements Box<T> {
  final Array<T> _entries;
  final (int, int) _strides;
  final (int, int) dimensions;
  final ({int x, int y}) offset;

  int get ld => _strides.$1;

  Matrix(
    int m,
    int n, {
    this.offset = oneIndexedMatrixOffset,
  })  : dimensions = (m, n),
        _strides = (m, 1),
        _entries = _Array<T>(m * n, offset: 0);

  Matrix.fromList(
    List<List<T>> list, {
    this.offset = oneIndexedMatrixOffset,
  })  : _entries = _Array<T>.fromList([
          for (var j = 0; j < (list.firstOrNull ?? []).length; j++) ...[
            for (var i = 0; i < list.length; i++) list[i][j],
          ],
        ], offset: 0),
        dimensions = (list.length, (list.firstOrNull ?? []).length),
        _strides = (list.length, 1);

  Matrix.fromData(
    List<T> data,
    this.dimensions, {
    this.offset = oneIndexedMatrixOffset,
  })  : assert(dimensions.$1 * dimensions.$2 >= data.length),
        _entries = Array.fromData(data, offset: 0),
        _strides = (dimensions.$1, 1);

  Matrix.fromSlice(
    Array<T> entries,
    this.dimensions, {
    this.offset = oneIndexedMatrixOffset,
  })  : //assert(dimensions.$1 * dimensions.$2 >= entries.length),
        _entries = _Array.fromData(entries.toData(), offset: 0),
        _strides = (dimensions.$1, 1);

  Matrix._(
    this._entries,
    this.dimensions,
    this._strides, {
    this.offset = oneIndexedMatrixOffset,
  });

  Matrix<T> call(
    int i,
    int j, {
    int? ld,
    ({int x, int y})? offset,
  }) {
    final entries = _entries(_getIndex(i, j), offset: 0);

    return Matrix._(
      entries,
      (ld ?? this.ld, dimensions.$2),
      (ld ?? this.ld, 1),
      offset: offset ?? this.offset,
    );
  }

  @pragma('vm:prefer-inline')
  int _getIndex(int i, int j) {
    i += this.offset.y;
    j += this.offset.x;
    return i * _strides.$2 + j * _strides.$1;
  }

  MatrixItemAccessor<T> operator [](int i) => MatrixItemAccessor(this, i);

  List<T> toData() => _entries.toData();

  Box<T> box(int i, int j) => this(i, j);

  Array<T> asArray() => Array.fromSlice(_entries.toData(), offset: offset.x);

  T get first => _entries.first;

  set first(T value) => _entries.first = value;

  Matrix<T> having({int? ld, int? lastd, ({int x, int y})? offset}) => Matrix._(
        _entries,
        (
          ld ?? this.dimensions.$1,
          ld != null && ld > 0 ? length ~/ ld : this.dimensions.$2
        ),
        (ld ?? this._strides.$1, this._strides.$2),
        offset: offset ?? this.offset,
      );

  Matrix<T> copy() => Matrix.fromData(toData(), dimensions, offset: offset);

  int get length => _entries.length;

  @override
  T get value => first;

  @override
  set value(T value) => first = value;

  Matrix<R> cast<R>() {
    final entries = _entries.cast<R>();
    return Matrix.fromSlice(entries, dimensions, offset: offset);
  }
}

class MatrixItemAccessor<T> {
  final int _i;
  final Matrix<T> _m;

  const MatrixItemAccessor(this._m, this._i);

  @pragma('vm:prefer-inline')
  T operator [](int j) => _m._entries[_m._getIndex(_i, j)];

  @pragma('vm:prefer-inline')
  void operator []=(int j, T value) => _m._entries[_m._getIndex(_i, j)] = value;
}

class _Array<T> implements Array<T> {
  @override
  final int offset;
  final List<T> _elements;

  const _Array.fromData(
    this._elements, {
    this.offset = oneIndexedArrayOffset,
  });

  @override
  _Array(int length, {this.offset = oneIndexedArrayOffset})
      : _elements = switch (T) {
          const (double) => Float64List(length),
          const (int) => Int64List(length),
          const (Complex) => Complex64List(length),
          const (bool) => List.filled(length, false),
          _ => throw UnimplementedError(),
        } as List<T>;

  _Array.fromList(List<T> list, {this.offset = oneIndexedArrayOffset})
      : _elements = switch (T) {
          const (double) => Float64List.fromList(list as List<double>),
          const (int) => Int64List.fromList(list as List<int>),
          const (Complex) => Complex64List.fromList(list as List<Complex>),
          const (bool) => [...list],
          _ => throw UnimplementedError(),
        } as List<T>;

  @override
  void assign(Array<T> array) {
    final dest = toData();
    final source = array.toData();
    dest.setRange(0, min(dest.length, source.length), source);
  }

  @override
  Array<T> slice(int index, {int? offset}) =>
      _slice(index, offset: offset ?? this.offset);

  _Array<T> _slice(int index, {int? offset, int? length}) {
    return _Array.fromData(
      switch (T) {
        const (double) => (_elements as Float64List).buffer.asFloat64List(
              (_elements as Float64List).offsetInBytes +
                  (index + this.offset) *
                      (_elements as Float64List).elementSizeInBytes,
              length,
            ),
        const (int) => (_elements as Int64List).buffer.asInt64List(
              (_elements as Int64List).offsetInBytes +
                  (index + this.offset) *
                      (_elements as Int64List).elementSizeInBytes,
              length,
            ),
        const (Complex) =>
          (_elements as Complex64List).slice(index + this.offset, length),
        const (bool) => _elements is ListSlice
            ? _elements.slice(
                index + this.offset,
                length != null
                    ? index + this.offset + length
                    : _elements.length,
              )
            : ListSlice(
                _elements,
                index + this.offset,
                length != null
                    ? index + this.offset + length
                    : _elements.length,
              ),
        _ => throw UnimplementedError(),
      } as List<T>,
      offset: offset ?? this.offset,
    );
  }

  @override
  @pragma('vm:prefer-inline')
  Array<T> call(int index, {int? offset}) {
    return slice(index, offset: offset);
  }

  @override
  @pragma('vm:prefer-inline')
  T operator [](int index) {
    return _elements[index + offset];
  }

  @override
  @pragma('vm:prefer-inline')
  void operator []=(int index, T value) {
    _elements[index + offset] = value;
  }

  @override
  Box<T> box(int index) {
    return DelegatingBox(() => this[index], (value) => this[index] = value);
  }

  @override
  List<T> toData() => _elements;

  @override
  Matrix<T> asMatrix([int ld = 0]) {
    return Matrix._(
      _Array.fromData(_elements, offset: 0),
      (ld, ld == 0 ? 0 : length ~/ ld),
      (ld, 1),
      offset: (x: offset, y: offset),
    );
  }

  @override
  Array<T> having({int? length, int? offset}) {
    return _slice(
      -this.offset,
      offset: offset ?? this.offset,
      length: length,
    );
  }

  @override
  Array<T> copy() => Array.fromData(toData(), offset: offset);

  @override
  Array<R> cast<R>() {
    if (T == R) return this as Array<R>;

    final typedData = switch (T) {
      const (bool) => throw UnimplementedError(),
      const (Complex) => (_elements as Complex64List).toData(),
      _ => _elements as TypedData,
    };

    final elements = switch (R) {
      const (double) => typedData.buffer.asFloat64List(typedData.offsetInBytes),
      const (int) => typedData.buffer.asInt64List(typedData.offsetInBytes),
      const (Complex) => Complex64List.fromData(
          typedData.buffer.asFloat64x2List(typedData.offsetInBytes)),
      _ => throw UnimplementedError(),
    } as List<R>;

    return _Array<R>.fromData(elements);
  }

  @override
  T get value => first;

  @override
  set value(T value) => first = value;

  @override
  int get length => _elements.length;

  @override
  T get first => _elements.first;

  @override
  set first(T value) => _elements.first = value;

  @override
  T get last => _elements.last;

  @override
  set last(T value) => _elements.last = value;

  @override
  List<T> operator +(List<T> other) =>
      _elements + (other is Array<T> ? other.toData() : other);

  @override
  void add(T value) => _elements.add(value);

  @override
  void addAll(Iterable<T> iterable) {
    iterable = iterable is Array<T> ? iterable.toData() : iterable;
    _elements.addAll(iterable);
  }

  @override
  bool any(bool Function(T element) test) => _elements.any(test);

  @override
  Map<int, T> asMap() => _elements.asMap();

  @override
  void clear() => _elements.clear();

  @override
  bool contains(Object? element) => _elements.contains(element);

  @override
  T elementAt(int index) => _elements.elementAt(index + offset);

  @override
  bool every(bool Function(T element) test) => _elements.every(test);

  @override
  Iterable<U> expand<U>(Iterable<U> Function(T element) toElements) =>
      _elements.expand(toElements);

  @override
  void fillRange(int start, int end, [T? fillValue]) =>
      _elements.fillRange(start + offset, end + offset, fillValue);

  @override
  T firstWhere(bool Function(T element) test, {T Function()? orElse}) =>
      _elements.firstWhere(test, orElse: orElse);

  @override
  U fold<U>(U initialValue, U Function(U previousValue, T element) combine) =>
      _elements.fold(initialValue, combine);

  @override
  Iterable<T> followedBy(Iterable<T> other) {
    other = other is Array<T> ? other.toData() : other;
    return _elements.followedBy(other);
  }

  @override
  void forEach(void Function(T element) action) => _elements.forEach(action);

  @override
  Iterable<T> getRange(int start, int end) =>
      _elements.getRange(start + offset, end + offset);

  @override
  int indexOf(T element, [int start = 0]) =>
      switch (_elements.indexOf(element, start + offset)) {
        -1 => -1,
        final i => i - offset,
      };

  @override
  int indexWhere(bool Function(T element) test, [int start = 0]) =>
      switch (_elements.indexWhere(test, start + offset)) {
        -1 => -1,
        final i => i - offset,
      };

  @override
  void insert(int index, T element) =>
      _elements.insert(index + offset, element);

  @override
  void insertAll(int index, Iterable<T> iterable) {
    iterable = iterable is Array<T> ? iterable.toData() : iterable;
    _elements.insertAll(index + offset, iterable);
  }

  @override
  bool get isEmpty => _elements.isEmpty;

  @override
  bool get isNotEmpty => _elements.isNotEmpty;

  @override
  Iterator<T> get iterator => _elements.iterator;

  @override
  String join([String separator = '']) => _elements.join(separator);

  @override
  int lastIndexOf(T element, [int? start]) => switch (_elements.lastIndexOf(
          element, start != null ? start + offset : null)) {
        -1 => -1,
        final i => i - offset,
      };

  @override
  int lastIndexWhere(bool Function(T element) test, [int? start]) =>
      switch (_elements.lastIndexWhere(
          test, start != null ? start + offset : null)) {
        -1 => -1,
        final i => i - offset,
      };

  @override
  T lastWhere(bool Function(T element) test, {T Function()? orElse}) =>
      _elements.lastWhere(test, orElse: orElse);

  @override
  set length(int newLength) => _elements.length;

  @override
  Iterable<U> map<U>(U Function(T e) toElement) => _elements.map(toElement);

  @override
  T reduce(T Function(T value, T element) combine) => _elements.reduce(combine);

  @override
  bool remove(Object? value) => _elements.remove(value);

  @override
  T removeAt(int index) => _elements.removeAt(index + offset);

  @override
  T removeLast() => _elements.removeLast();

  @override
  void removeRange(int start, int end) =>
      _elements.removeRange(start + offset, end + offset);

  @override
  void removeWhere(bool Function(T element) test) =>
      _elements.removeWhere((element) => false);

  @override
  void replaceRange(int start, int end, Iterable<T> replacements) {
    replacements =
        replacements is Array<T> ? replacements.toData() : replacements;
    _elements.replaceRange(start + offset, end + offset, replacements);
  }

  @override
  void retainWhere(bool Function(T element) test) =>
      _elements.retainWhere(test);

  @override
  Iterable<T> get reversed => _elements.reversed;

  @override
  void setAll(int index, Iterable<T> iterable) {
    iterable = iterable is Array<T> ? iterable.toData() : iterable;
    _elements.setAll(index + offset, iterable);
  }

  @override
  void setRange(int start, int end, Iterable<T> iterable, [int skipCount = 0]) {
    iterable = iterable is Array<T> ? iterable.toData() : iterable;
    _elements.setRange(start + offset, end + offset, iterable, skipCount);
  }

  @override
  void shuffle([Random? random]) => _elements.shuffle();

  @override
  T get single => _elements.single;

  @override
  T singleWhere(bool Function(T element) test, {T Function()? orElse}) =>
      _elements.singleWhere(test, orElse: orElse);

  @override
  Iterable<T> skip(int count) => _elements.skip(count);

  @override
  Iterable<T> skipWhile(bool Function(T value) test) =>
      _elements.skipWhile(test);

  @override
  void sort([int Function(T a, T b)? compare]) => _elements.sort();

  @override
  List<T> sublist(int start, [int? end]) =>
      _elements.sublist(start + offset, end != null ? end + offset : null);

  @override
  Iterable<T> take(int count) => _elements.take(count);

  @override
  Iterable<T> takeWhile(bool Function(T value) test) =>
      _elements.takeWhile(test);

  @override
  List<T> toList({bool growable = true}) =>
      _elements.toList(growable: growable);

  @override
  Set<T> toSet() => _elements.toSet();

  @override
  Iterable<T> where(bool Function(T element) test) => _elements.where(test);

  @override
  Iterable<U> whereType<U>() => _elements.whereType<U>();
}

class Matrix3d<T> {
  final Array<T> _entries;
  final (int, int, int) dimensions;
  final (int, int, int) _strides;
  final ({int x, int y, int z}) offset;

  Matrix3d(
    int m,
    int n,
    int o, {
    this.offset = oneIndexedMatrix3dOffset,
  })  : dimensions = (m, n, o),
        _strides = (m, n, 1),
        _entries = _Array<T>(m * n * o, offset: 0);

  Matrix3d.fromList(
    List<List<List<T>>> list, {
    this.offset = oneIndexedMatrix3dOffset,
  })  : _entries = _Array<T>.fromList(
          [
            for (var k = 0;
                k < (list.firstOrNull?.firstOrNull ?? []).length;
                k++) ...[
              for (var j = 0; j < (list.firstOrNull ?? []).length; j++) ...[
                for (var i = 0; i < list.length; i++) list[i][j][k],
              ],
            ],
          ],
          offset: 0,
        ),
        dimensions = (
          list.length,
          list[0].length,
          (list.firstOrNull?.firstOrNull ?? []).length,
        ),
        _strides = (list.length, list[0].length, 1);

  Matrix3d.fromData(
    List<T> data,
    this.dimensions, {
    this.offset = oneIndexedMatrix3dOffset,
  })  : assert(data.length >= dimensions.$1 * dimensions.$2 * dimensions.$3),
        _entries = _Array<T>.fromList(data, offset: 0),
        _strides = (dimensions.$1, dimensions.$2, 1);

  Matrix3d.fromSlice(
    Array<T> entries,
    this.dimensions, {
    this.offset = oneIndexedMatrix3dOffset,
  })  : //assert(entries.length >= dimensions.$1 * dimensions.$2 * dimensions.$3),
        _entries = _Array.fromData(entries.toData(), offset: 0),
        _strides = (dimensions.$1, dimensions.$2, 1);

  Matrix3d._(this._entries, this.dimensions, this._strides, this.offset);

  Matrix3d<T> call(
    int i,
    int j,
    int k, {
    ({int x, int y, int z}) offset = oneIndexedMatrix3dOffset,
  }) {
    final entries = _entries(_getIndex(i, j, k), offset: 0);
    return Matrix3d.fromSlice(entries, dimensions, offset: offset);
  }

  Matrix3dItemAccessor1<T> operator [](int i) => Matrix3dItemAccessor1(this, i);

  @pragma('vm:prefer-inline')
  int _getIndex(int i, int j, int k) {
    i += this.offset.y;
    j += this.offset.x;
    k += this.offset.z;
    final ld1 = _strides.$1;
    final ld2 = _strides.$2;
    return i + j * ld1 + k * ld2 * ld1;
  }

  List<T> toData() => _entries.toData();

  Box<T> box(int i, int j, int k) =>
      DelegatingBox(() => this[i][j][k], (value) => this[i][j][k] = value);

  Array<T> asArray() => _Array.fromData(_entries.toData(), offset: offset.x);

  Matrix3d<T> having({(int, int)? ld, ({int x, int y, int z})? offset}) =>
      Matrix3d._(
        _entries,
        (
          ld != null ? ld.$1 : this.dimensions.$1,
          ld != null ? ld.$2 : this.dimensions.$2,
          this.dimensions.$3,
        ),
        (
          ld != null ? ld.$1 : this._strides.$1,
          ld != null ? ld.$2 : this._strides.$2,
          this._strides.$3,
        ),
        offset ?? this.offset,
      );
}

class Matrix3dItemAccessor1<T> {
  final int _i;
  final Matrix3d<T> _m;

  const Matrix3dItemAccessor1(this._m, this._i);

  Matrix3dItemAccessor2<T> operator [](int j) =>
      Matrix3dItemAccessor2<T>(_m, _i, j);
}

class Matrix3dItemAccessor2<T> {
  final int _i;
  final int _j;
  final Matrix3d<T> _m;

  const Matrix3dItemAccessor2(this._m, this._i, this._j);

  @pragma('vm:prefer-inline')
  T operator [](int k) => _m._entries[_m._getIndex(_i, _j, k)];

  @pragma('vm:prefer-inline')
  void operator []=(int k, T value) =>
      _m._entries[_m._getIndex(_i, _j, k)] = value;
}
