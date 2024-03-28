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

  MatrixItemAcessor<T> operator [](int i) {
    return MatrixItemAcessor(this, i);
  }

  List<T> toData() {
    return _entries.toData();
  }

  Box<T> box(int i, int j) => this(i, j);

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

class MatrixItemAcessor<T> {
  final int _i;
  final Matrix<T> _m;

  const MatrixItemAcessor(this._m, this._i);

  T operator [](int j) => _m(_i, j).first;

  void operator []=(int j, T value) => _m(_i, j).first = value;
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
    return DelegatingBox(() => this[index], (value) => this[index] = value);
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
  void addAll(Iterable<T> iterable) => _elements.addAll(iterable);

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
      _elements.fillRange(start, end, fillValue);

  @override
  T firstWhere(bool Function(T element) test, {T Function()? orElse}) =>
      _elements.firstWhere(test, orElse: orElse);

  @override
  U fold<U>(U initialValue, U Function(U previousValue, T element) combine) =>
      _elements.fold(initialValue, combine);

  @override
  Iterable<T> followedBy(Iterable<T> other) => _elements.followedBy(other);

  @override
  void forEach(void Function(T element) action) => _elements.forEach(action);

  @override
  Iterable<T> getRange(int start, int end) => _elements.getRange(start, end);

  @override
  int indexOf(T element, [int start = 0]) =>
      _elements.indexOf(element, start + offset) - offset;

  @override
  int indexWhere(bool Function(T element) test, [int start = 0]) =>
      _elements.indexWhere(test, start + offset) - offset;

  @override
  void insert(int index, T element) =>
      _elements.insert(index + offset, element);

  @override
  void insertAll(int index, Iterable<T> iterable) =>
      _elements.insertAll(index + offset, iterable);

  @override
  bool get isEmpty => _elements.isEmpty;

  @override
  bool get isNotEmpty => _elements.isNotEmpty;

  @override
  Iterator<T> get iterator => _elements.iterator;

  @override
  String join([String separator = '']) => _elements.join(separator);

  @override
  int lastIndexOf(T element, [int? start]) =>
      _elements.lastIndexOf(element, start != null ? start + offset : null) -
      offset;

  @override
  int lastIndexWhere(bool Function(T element) test, [int? start]) =>
      _elements.lastIndexWhere(test, start != null ? start + offset : null) -
      offset;

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
  void replaceRange(int start, int end, Iterable<T> replacements) =>
      _elements.replaceRange(start + offset, end + offset, replacements);

  @override
  void retainWhere(bool Function(T element) test) =>
      _elements.retainWhere(test);

  @override
  Iterable<T> get reversed => _elements.reversed;

  @override
  void setAll(int index, Iterable<T> iterable) =>
      _elements.setAll(index + offset, iterable);

  @override
  void setRange(int start, int end, Iterable<T> iterable, [int skipCount = 0]) {
    final source = iterable is Array<T> ? iterable.toData() : iterable;
    _elements.setRange(start + offset, end + offset, source);
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
      _elements.takeWhile((value) => false);

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
