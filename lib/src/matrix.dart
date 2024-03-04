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

  factory Array.fromSlice(
    List<T> elements, {
    int offset = oneIndexedArrayOffset,
  }) {
    return _Array._(elements, offset: offset);
  }

  Array<T> slice(int index, {int? offset});

  Array<T> call(int index, {int? offset});

  T operator [](int index);
  void operator []=(int index, T value);

  Box<T> box(int index);

  List<T> toData();

  Matrix<T> asMatrix([int ld]);

  Array<T> dim([int? ld]);

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

typedef MatrixIndexer = int Function(List<int> dimentions, List<int> indexes);

int columnIndexed(List<int> dimentions, List<int> indexes) {
  dimentions = dimentions.sublist(0, indexes.length);
  var index = 0;
  var dim = 1;
  final dimIter = dimentions.reversed.iterator;
  for (final i in indexes) {
    dimIter.moveNext();
    dim *= dimIter.current;
    index += i * dim;
  }
  return index;
}

int rowIndexed(List<int> dimentions, List<int> indexes) {
  return columnIndexed(dimentions, indexes.reversed.toList());
}

const defaultMatrixIndexer = columnIndexed;

class Matrix<T> implements Box<T> {
  final Array<T> _entries;
  final (int, int) dimensions;
  final ({int x, int y}) offset;
  final MatrixIndexer _indexer;

  int get ld => dimensions.$1;

  Matrix(
    int m,
    int n, {
    this.offset = oneIndexedMatrixOffset,
    MatrixIndexer indexer = defaultMatrixIndexer,
  })  : dimensions = (m, 1),
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
        dimensions = (list.length, 1);

  Matrix.fromData(
    List<T> data,
    this.dimensions, {
    this.offset = oneIndexedMatrixOffset,
    MatrixIndexer indexer = defaultMatrixIndexer,
  })  : _entries = _Array.fromList(data, offset: 0),
        _indexer = indexer;

  Matrix.fromSlice(
    Array<T> entries,
    this.dimensions, {
    this.offset = oneIndexedMatrixOffset,
    MatrixIndexer indexer = defaultMatrixIndexer,
  })  : _entries = _Array._(entries.toData(), offset: 0),
        _indexer = indexer;

  Matrix<T> call(
    int i,
    int j, {
    int? ld,
    ({int x, int y})? offset,
  }) {
    var entries = _entries(
      _indexer([this.dimensions.$1, this.dimensions.$2],
          [i + this.offset.y, j + this.offset.x]),
      offset: 0,
    );

    return Matrix.fromSlice(
      entries,
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

  Matrix<T> dim(int ld) =>
      Matrix.fromSlice(_entries, (ld, 1), offset: offset, indexer: _indexer);

  @override
  T get value => first;

  @override
  set value(T value) => first = value;
}

class _MatrixArrayAdapter<T> implements Array<T> {
  final int i;
  final Matrix<T> _m;

  const _MatrixArrayAdapter(this._m, this.i);

  @override
  Array<T> slice(int j, {int? offset}) {
    final slice = _m(i, j);
    return Array.fromSlice(slice._entries.toData(),
        offset: offset ?? _m.offset.x);
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
    return Matrix.fromSlice(
      getEntries(),
      ld != null ? (ld, 1) : _m.dimensions,
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
  Array<T> dim([int? ld]) => getEntries().dim(ld);

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
  final int offset;
  final List<T> _elements;

  _Array._(
    this._elements, {
    this.offset = oneIndexedArrayOffset,
  });

  @override
  _Array(int length, {this.offset = oneIndexedArrayOffset})
      : _elements = switch (T) {
          double => Float64List(length),
          int => Int64List(length),
          Complex => Float64x2List(length),
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
  Array<T> slice(int index, {int? offset}) =>
      _slice(index, offset: offset ?? this.offset);

  _Array<T> _slice(int index, {int? offset}) {
    return _Array._(
      switch (T) {
        double => (_elements as Float64List).buffer.asFloat64List(
              (_elements as Float64List).offsetInBytes +
                  (index + this.offset) *
                      (_elements as Float64List).elementSizeInBytes,
            ),
        int => (_elements as Int64List).buffer.asInt64List(
              (_elements as Int64List).offsetInBytes +
                  (index + this.offset) *
                      (_elements as Int64List).elementSizeInBytes,
            ),
        Complex => (_elements as Complex64List).slice(index + this.offset),
        bool => _elements is ListSlice
            ? _elements.slice(index + this.offset, _elements.length)
            : ListSlice(_elements, index + this.offset, _elements.length),
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
    return Matrix.fromSlice(
      _Array._(_elements, offset: 0),
      (ld, 1),
      offset: (x: offset, y: offset),
    );
  }

  @override
  Array<T> dim([int? ld]) => _slice(-offset, offset: offset);

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
  final ({int x, int y, int z}) offset;

  Matrix3d(
    int m,
    int n,
    int o, {
    this.offset = oneIndexedMatrix3dOffset,
  })  : dimensions = (m, n, 1),
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
        dimensions = (list.length, list[0].length, 1);

  Matrix3d.fromData(
    List<T> data,
    this.dimensions, {
    this.offset = oneIndexedMatrix3dOffset,
  }) : _entries = _Array<T>.fromList(data, offset: 0);

  Matrix3d.fromSlice(
    Array<T> entries,
    this.dimensions, {
    this.offset = oneIndexedMatrix3dOffset,
  }) : _entries = _Array._(entries.toData(), offset: 0);

  Matrix3d<T> call(
    int i,
    int j,
    int k, {
    ({int x, int y, int z}) offset = oneIndexedMatrix3dOffset,
  }) {
    final ld1 = dimensions.$1;
    final ld2 = dimensions.$2;
    final entries = _entries(
        (i + this.offset.x) +
            (j + this.offset.y) * ld1 +
            (k + this.offset.z) * ld2 * ld1,
        offset: 0);
    return Matrix3d.fromSlice(entries, dimensions, offset: offset);
  }

  Matrix<T> operator [](int i) {
    final ld1 = dimensions.$1;
    final ld2 = dimensions.$2;
    return Matrix.fromSlice(
      _entries(i + offset.x, offset: 0),
      (ld2, ld1),
      offset: (x: offset.y, y: offset.z),
    );
  }

  List<T> toData() {
    return _entries.toData();
  }

  Box<T> box(int i, int j, int k) => this[i][j].box(k);

  Array<T> asArray() => _Array._(_entries.toData(), offset: offset.x);
}
