import 'dart:math';
import 'dart:typed_data';

import 'package:collection/collection.dart';
import 'package:lapack/src/box.dart';
import 'package:lapack/src/complex.dart';

const oneIndexedArrayOffset = -1;
const zeroIndexedArrayOffset = 0;
const oneIndexedMatrixOffset = (x: -1, y: -1);
const zeroIndexedMatrixOffset = (x: 0, y: 0);
const oneIndexedMatrix3dOffset = (x: -1, y: -1, z: -1);
const zeroIndexedMatrix3dOffset = (x: 0, y: 0, z: 0);

// TODO: Make it iterable
// TODO: List interface
abstract interface class Array<T> implements Box<T> {
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

  T operator [](int index);
  void operator []=(int index, T value);

  Box<T> box(int index);

  List<T> toData();

  Matrix<T> asMatrix([int ld]);

  Array<T> having({int? length, int? offset});

  Array<T> copy();

  Array<R> cast<R>();

  int get offset;

  int get length;

  T get first;

  set first(T value);
}

extension ArrayExtension<T> on Array<T> {
  T maxval(int start, int end) {
    switch (T) {
      case double:
      case int:
        var value = this[start] as num;
        for (var i = start + 1; i <= end; i++) {
          if (value < (this[i] as num)) value = this[i] as num;
        }
        return value as T;

      case Complex:
        var value = this[start] as Complex;
        for (var i = start + 1; i <= end; i++) {
          if (value < (this[i] as Complex)) value = this[i] as Complex;
        }
        return value as T;

      case bool:
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
      case double:
      case int:
        int loc = start;
        var value = this[start] as num;
        for (var i = start + 1; i <= end; i++) {
          if (value < (this[i] as num)) {
            value = this[i] as num;
            loc = i;
          }
        }
        return loc;

      case Complex:
        int loc = start;
        Complex value = this[start] as Complex;
        for (var i = start + 1; i <= end; i++) {
          if (value < (this[i] as Complex)) {
            value = this[i] as Complex;
            loc = i;
          }
        }
        return loc;

      case bool:
        for (var i = start; i <= end; i++) {
          if (this[start] as bool) return i;
        }
        return 0;
      default:
        throw UnimplementedError();
    }
  }
}

typedef MatrixIndexer = int Function(List<int> strides, List<int> indexes);

int columnIndexed(List<int> strides, List<int> indexes) {
  strides = strides.sublist(0, indexes.length);
  var index = 0;
  var stride = 1;
  final strideIter = strides.reversed.iterator;
  for (final i in indexes) {
    strideIter.moveNext();
    stride *= strideIter.current;
    index += i * stride;
  }
  return index;
}

int rowIndexed(List<int> strides, List<int> indexes) {
  return columnIndexed(strides, indexes.reversed.toList());
}

const defaultMatrixIndexer = columnIndexed;

class Matrix<T> implements Box<T> {
  final Array<T> _entries;
  final (int, int) _strides;
  final (int, int) dimensions;
  final ({int x, int y}) offset;
  final MatrixIndexer _indexer;

  int get ld => _strides.$1;

  Matrix(
    int m,
    int n, {
    this.offset = oneIndexedMatrixOffset,
    MatrixIndexer indexer = defaultMatrixIndexer,
  })  : dimensions = (m, n),
        _strides = (m, 1),
        _indexer = indexer,
        _entries = _Array<T>(m * n, offset: 0);

  Matrix.fromList(
    List<List<T>> list, {
    this.offset = oneIndexedMatrixOffset,
    MatrixIndexer indexer = defaultMatrixIndexer,
  })  : _entries = _Array<T>.fromList([
          for (var j = 0; j < (list.firstOrNull ?? []).length; j++) ...[
            for (var i = 0; i < list.length; i++) list[i][j],
          ],
        ], offset: 0),
        _indexer = indexer,
        dimensions = (list.length, (list.firstOrNull ?? []).length),
        _strides = (list.length, 1);

  Matrix.fromData(
    List<T> data,
    this.dimensions, {
    this.offset = oneIndexedMatrixOffset,
    MatrixIndexer indexer = defaultMatrixIndexer,
  })  : assert(dimensions.$1 * dimensions.$2 >= data.length),
        _entries = Array.fromData(data, offset: 0),
        _indexer = indexer,
        _strides = (dimensions.$1, 1);

  Matrix.fromSlice(
    Array<T> entries,
    this.dimensions, {
    this.offset = oneIndexedMatrixOffset,
    MatrixIndexer indexer = defaultMatrixIndexer,
  })  : //assert(dimensions.$1 * dimensions.$2 >= entries.length),
        _entries = _Array.fromData(entries.toData(), offset: 0),
        _strides = (dimensions.$1, 1),
        _indexer = indexer;

  Matrix._(
    this._entries,
    this.dimensions,
    this._strides, {
    this.offset = oneIndexedMatrixOffset,
    MatrixIndexer indexer = defaultMatrixIndexer,
  }) : _indexer = indexer;

  Matrix<T> call(
    int i,
    int j, {
    int? ld,
    ({int x, int y})? offset,
  }) {
    var entries = _entries(
      _indexer([this._strides.$1, this._strides.$2],
          [i + this.offset.y, j + this.offset.x]),
      offset: 0,
    );

    return Matrix._(
      entries,
      (ld ?? this.ld, dimensions.$2),
      (ld ?? this.ld, 1),
      offset: offset ?? this.offset,
      indexer: _indexer,
    );
  }

  Array<T> operator [](int i) {
    return _MatrixArrayAdapter(this, i);
  }

  List<T> toData() {
    return _entries.toData();
  }

  Box<T> box(int i, int j) => this[i].box(j);

  Array<T> asArray() => Array.fromSlice(_entries.toData());

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
        indexer: _indexer,
      );

  Matrix<T> copy() {
    return Matrix.fromData(toData(), dimensions,
        offset: offset, indexer: _indexer);
  }

  int get length => _entries.length;

  @override
  T get value => first;

  @override
  set value(T value) => first = value;

  Matrix<R> cast<R>() {
    final entries = _entries.cast<R>();
    return Matrix.fromSlice(entries, dimensions,
        offset: offset, indexer: _indexer);
  }
}

class _MatrixArrayAdapter<T> implements Array<T> {
  final int i;
  final Matrix<T> _m;

  const _MatrixArrayAdapter(this._m, this.i);

  @override
  void assign(Array<T> array) {
    final entries = getEntries();
    final dest = entries.toData();
    final source = array.toData();
    dest.setRange(0, min(dest.length, source.length), source);
  }

  @override
  Array<T> slice(int j, {int? offset}) {
    final slice = _m(i, j);
    return Array.fromSlice(slice._entries.toData(),
        offset: offset ?? this.offset);
  }

  @override
  Array<T> call(int j, {int? offset}) {
    return slice(j, offset: offset);
  }

  @override
  Box<T> box(int j) {
    return slice(j);
  }

  Array<T> getEntries([int? j]) {
    return _m(i, j ?? -_m.offset.x)._entries;
  }

  @override
  List<T> toData() {
    return getEntries().toData();
  }

  @override
  Matrix<T> asMatrix([int? ld]) {
    return Matrix._(
      _Array.fromData(getEntries().toData(), offset: 0),
      ld != null ? (ld, _m.dimensions.$2) : _m.dimensions,
      ld != null ? (ld, 1) : _m._strides,
      offset: _m.offset,
    );
  }

  @override
  T operator [](int j) {
    return _m(i, j).first;
  }

  @override
  void operator []=(int j, T value) {
    _m(i, j).first = value;
  }

  @override
  Array<T> having({int? length, int? offset}) {
    final entries = getEntries();
    return entries.having(
      length: length,
      offset: offset ?? this.offset,
    );
  }

  @override
  Array<T> copy() {
    final entries = getEntries();
    return entries.copy();
  }

  @override
  Array<R> cast<R>() {
    return getEntries().cast<R>();
  }

  @override
  int get offset => _m.offset.x;

  @override
  int get length => getEntries().length;

  @override
  T get value => first;

  @override
  set value(T value) => first = value;

  @override
  T get first => getEntries().first;

  @override
  set first(T value) => getEntries().first = value;
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
          double => Float64List(length),
          int => Int64List(length),
          Complex => Complex64List(length),
          bool => List.filled(length, false),
          _ => throw UnimplementedError(),
        } as List<T>;

  _Array.fromList(List<T> list, {this.offset = oneIndexedArrayOffset})
      : _elements = switch (T) {
          double => Float64List.fromList(list as List<double>),
          int => Int64List.fromList(list as List<int>),
          Complex => Complex64List.fromList(list as List<Complex>),
          bool => [...list],
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
        double => (_elements as Float64List).buffer.asFloat64List(
              (_elements as Float64List).offsetInBytes +
                  (index + this.offset) *
                      (_elements as Float64List).elementSizeInBytes,
              length,
            ),
        int => (_elements as Int64List).buffer.asInt64List(
              (_elements as Int64List).offsetInBytes +
                  (index + this.offset) *
                      (_elements as Int64List).elementSizeInBytes,
              length,
            ),
        Complex =>
          (_elements as Complex64List).slice(index + this.offset, length),
        bool => _elements is ListSlice
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
  Array<T> call(int index, {int? offset}) {
    return slice(index, offset: offset);
  }

  @override
  T operator [](int index) {
    return _elements[index + offset];
  }

  @override
  void operator []=(int index, T value) {
    _elements[index + offset] = value;
  }

  @override
  Box<T> box(int index) {
    return _ArrayElementBox(this, index);
  }

  @override
  List<T> toData() {
    return _elements;
  }

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
  Array<T> copy() {
    return Array.fromData(toData(), offset: offset);
  }

  @override
  Array<R> cast<R>() {
    if (T == R) return this as Array<R>;

    final typedData = switch (T) {
      bool => throw UnimplementedError(),
      Complex => (_elements as Complex64List).toData(),
      _ => _elements as TypedData,
    };

    final elements = switch (R) {
      double => typedData.buffer.asFloat64List(typedData.offsetInBytes),
      int => typedData.buffer.asInt64List(typedData.offsetInBytes),
      Complex => Complex64List.fromData(
          typedData.buffer.asFloat64x2List(typedData.offsetInBytes)),
      _ => throw UnimplementedError(),
    } as List<R>;

    return _Array<R>.fromData(elements);
  }

  @override
  int get length => _elements.length;

  @override
  T get value => first;

  @override
  set value(T value) => first = value;

  @override
  T get first => _elements.first;

  @override
  set first(T value) => _elements.first = value;
}

class _ArrayElementBox<T> implements Box<T> {
  final Array<T> _array;
  final int _index;
  _ArrayElementBox(this._array, this._index);

  @override
  T get value => _array[_index];
  @override
  set value(T value) => _array[_index] = value;

  @override
  String toString() => value.toString();
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

  Matrix3d._(
    this._entries,
    this.dimensions,
    this._strides,
    this.offset,
  );

  Matrix3d<T> call(
    int i,
    int j,
    int k, {
    ({int x, int y, int z}) offset = oneIndexedMatrix3dOffset,
  }) {
    final ld1 = _strides.$1;
    final ld2 = _strides.$2;
    final entries = _entries(
        (i + this.offset.x) +
            (j + this.offset.y) * ld1 +
            (k + this.offset.z) * ld2 * ld1,
        offset: 0);
    return Matrix3d.fromSlice(entries, dimensions, offset: offset);
  }

  Matrix<T> operator [](int i) {
    final ld1 = _strides.$1;
    final ld2 = _strides.$2;
    return Matrix._(
      _entries(i + offset.x, offset: 0),
      (ld2, ld1),
      (ld2, ld1),
      offset: (x: offset.y, y: offset.z),
    );
  }

  List<T> toData() {
    return _entries.toData();
  }

  Box<T> box(int i, int j, int k) => this[i][j].box(k);

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
