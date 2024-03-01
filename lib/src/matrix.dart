import 'dart:typed_data';

import 'package:collection/collection.dart';
import 'package:lapack/src/box.dart';

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
    return _Array.fromSlice(elements, offset: offset);
  }

  Array<T> slice(
    int index, {
    int offset = oneIndexedArrayOffset,
  });

  Array<T> call(
    int index, {
    int offset = oneIndexedArrayOffset,
  });

  T operator [](int index);
  void operator []=(int index, T value);

  Box<T> box(int index);

  List<T> toRawList();

  Matrix<T> asMatrix([int ld]);

  T maxval(int start, int end);
  int maxloc(int start, int end, {int dim = 1});

  Array<T> dim([int? ld]);

  T get first;

  set first(T value);
}

typedef MatrixIndexer = int Function(int ld, int i, int j);

int columnIndexed(int ld, int i, int j) {
  return j * ld + i;
}

int rowIndexed(int ld, int i, int j) {
  return j + i * ld;
}

const defaultMatrixIndexer = columnIndexed;

class Matrix<T> implements Box<T> {
  final Array<T> _entries;
  final int _ld;
  final ({int x, int y}) offset;
  final MatrixIndexer _indexer;

  Matrix(
    int m,
    int n, {
    this.offset = oneIndexedMatrixOffset,
    MatrixIndexer indexer = defaultMatrixIndexer,
  })  : _ld = m,
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
        _ld = list.length;

  Matrix.fromSlice(
    Array<T> entries,
    this._ld, {
    this.offset = oneIndexedMatrixOffset,
    MatrixIndexer indexer = defaultMatrixIndexer,
  })  : _entries = _Array.fromSlice(entries.toRawList(), offset: 0),
        _indexer = indexer;

  Matrix<T> call(int i, int j, {int? ld}) {
    var entries = _entries(_indexer(this._ld, i + offset.y, j + offset.x));
    return Matrix.fromSlice(
      entries,
      ld ?? _ld,
      indexer: _indexer,
    );
  }

  Array<T> operator [](int i) {
    return _MatrixArrayAdapter(this, i);
  }

  List<T> toRawList() {
    return _entries.toRawList();
  }

  Box<T> box(int i, int j) => this[i].box(j);

  Array<T> asArray() => Array.fromSlice(_entries.toRawList());

  T get first => _entries.first;

  set first(T value) => _entries.first = value;

  Matrix<T> dim(int ld) => this(-offset.y, -offset.x, ld: ld);

  @override
  T get value => first;

  @override
  set value(T value) => first = value;
}

class _MatrixArrayAdapter<T> implements Array<T> {
  final Array<T> _entries;
  final Matrix<T> _m;

  _MatrixArrayAdapter(this._m, int i)
      : _entries = _Array.fromSlice(
          _m._entries
              .slice(_m._indexer(_m._ld, i + _m.offset.y, 0))
              .toRawList(),
          offset: 0,
        );

  @override
  Array<T> slice(
    int j, {
    int offset = oneIndexedArrayOffset,
  }) {
    return _entries(_getIndex(j), offset: offset);
  }

  @override
  Array<T> call(
    int j, {
    int offset = oneIndexedArrayOffset,
  }) {
    return slice(j, offset: offset);
  }

  @override
  Box<T> box(int j) {
    return slice(j);
  }

  @override
  List<T> toRawList() {
    return _entries.toRawList();
  }

  @override
  Matrix<T> asMatrix([int? ld]) {
    return Matrix.fromSlice(
      _entries,
      ld ?? _m._ld,
      offset: _m.offset,
    );
  }

  @override
  T maxval(int start, int end) {
    return _entries.maxval(start, end);
  }

  @override
  int maxloc(int start, int end, {int dim = 1}) {
    return _entries.maxloc(start, end, dim: dim);
  }

  int _getIndex(int j) {
    j += _m.offset.x;
    final i = j ~/ _m._ld;
    j -= i * _m._ld;
    return _m._indexer(_m._ld, i, j);
  }

  @override
  T operator [](int j) {
    return _entries[_getIndex(j)];
  }

  @override
  void operator []=(int j, T value) {
    _entries[_getIndex(j)] = value;
  }

  @override
  Array<T> dim([int? ld]) => this;

  @override
  T get value => first;

  @override
  set value(T value) => first = value;

  @override
  T get first => _entries.first;

  @override
  set first(T value) => _entries.first = value;
}

class _Array<T> implements Array<T> {
  final int offset;
  final List<T> _elements;

  @override
  _Array(
    int length, {
    this.offset = oneIndexedArrayOffset,
  }) : _elements = switch (T) {
          double => Float64List(length),
          int => Int64List(length),
          bool => List.filled(length, false),
          _ => throw UnimplementedError(),
        } as List<T>;

  _Array.fromList(
    List<T> list, {
    this.offset = oneIndexedArrayOffset,
  }) : _elements = switch (T) {
          double => Float64List.fromList(list as List<double>),
          int => Int64List.fromList(list as List<int>),
          bool => [...list],
          _ => throw UnimplementedError(),
        } as List<T>;

  _Array.fromSlice(
    this._elements, {
    this.offset = oneIndexedArrayOffset,
  });

  @override
  Array<T> slice(
    int index, {
    int offset = oneIndexedArrayOffset,
  }) {
    return _Array.fromSlice(
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
        bool => _elements is ListSlice
            ? _elements.slice(index + this.offset, _elements.length)
            : ListSlice(_elements, index + this.offset, _elements.length),
        _ => throw UnimplementedError(),
      } as List<T>,
      offset: offset,
    );
  }

  @override
  Array<T> call(
    int index, {
    int offset = oneIndexedArrayOffset,
  }) {
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
  List<T> toRawList() {
    return _elements;
  }

  @override
  Matrix<T> asMatrix([int ld = 0]) {
    return Matrix.fromSlice(
      _Array.fromSlice(_elements, offset: 0),
      ld,
      offset: (x: offset, y: offset),
    );
  }

  @override
  T maxval(int start, int end) {
    switch (T) {
      case double:
      case int:
        var value = this[start] as num;
        for (var i = start + 1; i <= end; i++) {
          if (value < (this[i] as num)) value = this[i] as num;
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

  @override
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

      case bool:
        for (var i = start; i <= end; i++) {
          if (this[start] as bool) return i;
        }
        return 0;
      default:
        throw UnimplementedError();
    }
  }

  @override
  Array<T> dim([int? ld]) => this;

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
  final (int, int) _ld;
  final ({int x, int y, int z}) offset;

  Matrix3d(
    int m,
    int n,
    int o, {
    this.offset = oneIndexedMatrix3dOffset,
  })  : _ld = (m, n),
        _entries = _Array<T>(m * n * o, offset: 0);

  Matrix3d.fromList(
    List<List<List<T>>> list, {
    this.offset = oneIndexedMatrix3dOffset,
  })  : _entries = _Array<T>.fromList(
          [
            for (var k = 0; k < (list.firstOrNull ?? []).length; k++) ...[
              for (var i = 0; i < (list.firstOrNull ?? []).length; i++) ...[
                for (var j = 0; j < list.length; j++) list[i][j][k],
              ],
            ],
          ],
          offset: 0,
        ),
        _ld = (list.length, list[0].length);

  Matrix3d.fromSlice(
    Array<T> entries,
    this._ld, {
    this.offset = oneIndexedMatrix3dOffset,
  }) : _entries = _Array.fromSlice(entries.toRawList(), offset: 0);

  Matrix3d<T> call(
    int i,
    int j,
    int k, {
    ({int x, int y, int z}) offset = oneIndexedMatrix3dOffset,
  }) {
    final entries = _entries((i + this.offset.x) * _ld.$1 + j + this.offset.y,
        offset: 0)(k + this.offset.z);
    return Matrix3d.fromSlice(entries, _ld, offset: offset);
  }

  Matrix<T> operator [](int i) {
    return Matrix.fromSlice(
      _entries((i + offset.x) * _ld.$1),
      _ld.$1 * _ld.$2,
      offset: (x: offset.y, y: offset.z),
    );
  }

  List<T> toRawList() {
    return _entries.toRawList();
  }

  Box<T> box(int i, int j, int k) => this[i][j].box(k);

  Array<T> asArray() => _Array.fromSlice(_entries.toRawList());
}
