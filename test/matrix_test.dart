// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/complex.dart';
import 'package:dart_lapack/src/matrix.dart';
import 'package:test/test.dart';

void main() {
  group('Array', () {
    test('fromList', () {
      final a = Array.fromList([1, 2, 3, 4]);
      expect(a[1], 1);
      expect(a[2], 2);
      expect(a[3], 3);
      expect(a[4], 4);
      expect(a.toData(), [1, 2, 3, 4]);
    });

    test('assign', () {
      final a = Array.fromList([1, 2, 3, 4]);
      final b = Array.fromList([0, 0, 0, 0]);

      // equal lengths
      b.assign(a);
      expect(b.toData(), [1, 2, 3, 4]);

      // dest smaller
      final c = b(3);
      c.assign(Array.fromList([9, 8, 7, 6]));
      expect(c.toData(), [9, 8]);
      expect(b.toData(), [1, 2, 9, 8]);

      // source smaller
      final d = Array.fromList([99, 88]);
      b.assign(d);
      expect(b.toData(), [99, 88, 9, 8]);
    });

    test('copy', () {
      final a = Array.fromList([1, 2, 3, 4]);
      final b = a.copy();
      expect(b.toData(), [1, 2, 3, 4]);
      b[1] = 999;
      expect(b.toData(), [999, 2, 3, 4]);
      expect(a.toData(), [1, 2, 3, 4]);
      final c = a(3);
      expect(c.toData(), [3, 4]);
      final d = c.copy();
      expect(d.toData(), [3, 4]);
    });

    test('offset', () {
      final a = Array.fromList([1, 2, 3, 4], offset: -10);
      expect(a[10], 1);
      expect(a[11], 2);
      expect(a[12], 3);
      expect(a[13], 4);
      expect(a.toData(), [1, 2, 3, 4]);
    });

    test('sliced offset', () {
      final a = Array.fromList([1, 2, 3, 4], offset: -10);
      final b = a.slice(12, offset: 0);
      expect(b[0], 3);
      expect(b[1], 4);
      expect(b.toData(), [3, 4]);
    });

    test('maxval', () {
      final a = Array.fromList([9, 1, 2, 3, 4]);
      expect(a.maxval(2, 5), 4);
      expect(a.maxval(2, 4), 3);
      expect(a.maxval(1, 4), 9);
    });

    test('maxloc', () {
      final a = Array.fromList([9, 1, 2, 3, 4]);
      expect(a.maxloc(2, 5), 4);
      expect(a.maxloc(2, 4), 3);
      expect(a.maxloc(1, 4), 1);
    });

    test('assignment', () {
      final a = Array.fromList([1, 2, 3, 4]);
      expect(a[2], 2);
      a[2] = 999;
      expect(a[1], 1);
      expect(a[2], 999);
      expect(a[3], 3);
      expect(a[4], 4);
    });

    test('slice', () {
      final a = Array.fromList([1, 2, 3, 4])(2);
      expect(a[1], 2);
      expect(a[2], 3);
      expect(a[3], 4);
      expect(a.toData(), [2, 3, 4]);

      final b = a(2);
      expect(b[1], 3);
      expect(b[2], 4);
      expect(b.toData(), [3, 4]);
    });

    test('sliced assignment', () {
      final a = Array.fromList([1, 2, 3, 4])(2);
      expect(a[2], 3);
      a[2] = 999;
      expect(a[1], 2);
      expect(a[2], 999);
      expect(a[3], 4);

      final b = a(2);
      expect(b[1], 999);
      a[2] = 123;
      expect(b[1], 123);
      expect(b[2], 4);
    });

    test('box', () {
      final a = Array.fromList([1, 2, 3, 4]);
      expect(a[3], 3);
      final box = a.box(3);
      box.value = 999;
      expect(a[1], 1);
      expect(a[2], 2);
      expect(a[3], 999);
      expect(a[4], 4);
    });

    test('sliced box', () {
      final a = Array.fromList([1, 2, 3, 4]).slice(2);
      expect(a[2], 3);
      final box = a.box(2);
      box.value = 999;
      expect(a[1], 2);
      expect(a[2], 999);
      expect(a[3], 4);
    });

    test('asMatrix', () {
      final a = Array.fromList([1, 2, 3, 4, 5, 6, 7, 8, 9]);
      final m = a.asMatrix(3);
      expect(m[1][1], 1);
      expect(m[1][2], 4);
      expect(m[1][3], 7);
      expect(m[2][1], 2);
      expect(m[2][2], 5);
      expect(m[2][3], 8);
      expect(m[3][1], 3);
      expect(m[3][2], 6);
      expect(m[3][3], 9);
    });

    test('sliced asMatrix', () {
      final a = Array.fromList([1, 2, 3, 4, 5, 6, 7, 8, 9]);
      final m = a(2).asMatrix(3);
      expect(m[1][1], 2);
      expect(m[1][2], 5);
      expect(m[1][3], 8);
      expect(m[2][1], 3);
      expect(m[2][2], 6);
      expect(m[2][3], 9);
      expect(m[3][1], 4);
      expect(m[3][2], 7);
      // expect(m[3][3], 0);
    });

    test('zero indexed', () {
      final a = Array.fromList([1, 2, 3, 4], offset: 0);
      expect(a[0], 1);
      expect(a[1], 2);
      expect(a[2], 3);
      expect(a[3], 4);
    });

    test('slice zero indexed', () {
      final a = Array.fromList([1, 2, 3, 4], offset: 0);
      expect(a[0], 1);
      expect(a[1], 2);
      expect(a[2], 3);
      expect(a[3], 4);

      final b = a(2);
      expect(b[0], 3);
      expect(b[1], 4);

      // back to one-indexed
      final c = a(2, offset: oneIndexedArrayOffset);
      expect(c[1], 3);
      expect(c[2], 4);
    });

    test('first', () {
      final a = Array.fromList([1, 2, 3, 4]);
      expect(a.first, 1);
      a.first = 999;
      expect(a[1], 999);
      expect(a[2], 2);
      expect(a[3], 3);
      expect(a[4], 4);
    });

    test('Complex array', () {
      final a = Array.fromList(
          [(1.0, 2.0), (2.0, 3.0), (3.0, 4.0), (4.0, 5.0)].toComplexList());
      expect(a.first, Complex(1.0, 2.0));
      expect(a[1], Complex(1, 2));
      expect(a[2], Complex(2, 3));
      expect(a[3], Complex(3, 4));
      expect(a[4], Complex(4, 5));
    });

    test('having', () {
      final a = Array.fromList([1, 2, 3, 4]);
      final b = a.having(offset: -10);
      expect([b[10], b[11], b[12], b[13]], [1, 2, 3, 4]);
      final c = b.having(length: 2);
      expect(a.length, 4);
      expect(c.length, 2);
      expect(c.first, 1);
      expect(c[10], 1);
    });

    test('cast', () {
      final a = Array.fromList([1.0, 2.0, 3.0, 4.0]);
      final b = a.cast<Complex>();
      expect(b.toData(), [Complex(1, 2), Complex(3, 4)]);
    });

    group('List', () {
      test('length', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.length, 4);

        final b = Array.fromList(<int>[]);
        expect(b.length, 0);

        final c = a(3);
        expect(c.length, 2);
      });

      test('first/last', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.first, 1);
        expect(a.last, 4);

        final b = a(3);
        expect(b.first, 3);
        expect(b.last, 4);
      });

      test('concatenation', () {
        final a = Array.fromList([1, 2, 3, 4]);
        final b = Array.fromList([5, 6, 7, 8]);

        final c = a + b;
        expect(c, [1, 2, 3, 4, 5, 6, 7, 8]);
      });

      test('add', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(() => a.add(5), throwsUnsupportedError);
      });

      test('addAll', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(() => a.addAll([5, 6, 7, 8]), throwsUnsupportedError);
      });

      test('any', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.any((e) => e == 3), true);
        expect(a.any((e) => e == 5), false);
      });

      test('asMap', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.asMap(), {0: 1, 1: 2, 2: 3, 3: 4});
      });

      test('clear', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(() => a.clear(), throwsUnsupportedError);
      });

      test('contains', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.contains(3), true);
        expect(a.contains(5), false);
      });

      test('elementAt', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.elementAt(2), 2);
        expect(() => a.elementAt(5), throwsRangeError);
      });

      test('every', () {
        final a = Array.fromList([2, 2, 2, 2]);
        expect(a.every((e) => e == 2), true);

        final b = Array.fromList([2, 2, 3, 2]);
        expect(b.every((e) => e == 2), false);
      });

      test('expand', () {
        Iterable<int> duplicate(int n) sync* {
          yield n;
          yield n;
        }

        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.expand(duplicate), [1, 1, 2, 2, 3, 3, 4, 4]);
      });

      test('fillRange', () {
        final a = Array.fromList([1, 2, 3, 4]);
        a.fillRange(2, 4, 99);
        expect(a, [1, 99, 99, 4]);
      });

      test('firstWhere', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.firstWhere((e) => e % 2 == 0), 2);
        expect(a.firstWhere((e) => e % 5 == 0, orElse: () => 99), 99);
      });

      test('fold', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.fold('', (r, n) => '$r$n'), '1234');

        final b = a(2);
        expect(b.fold('', (r, n) => '$r$n'), '234');

        final c = b.having(length: 2);
        expect(c.fold('', (r, n) => '$r$n'), '23');
      });

      test('followedBy', () {
        final a = Array.fromList([1, 2, 3, 4]);
        final b = Array.fromList([5, 6, 7, 8]);

        final c = a.followedBy(b);
        expect(c, [1, 2, 3, 4, 5, 6, 7, 8]);
      });

      test('forEach', () {
        final a = Array.fromList([1, 2, 3, 4]);

        final r = <int>[];
        void push(int n) => r.add(n);
        a.forEach(push);
        expect(r, [1, 2, 3, 4]);
      });

      test('getRange', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.getRange(2, 4), [2, 3]);
      });

      test('indexOf', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.indexOf(2), 2);
        expect(a.indexOf(99), -1);
        expect(a.indexOf(2, 2), 2);
        expect(a.indexOf(2, 3), -1);
      });

      test('indexWhere', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.indexWhere((e) => e == 2), 2);
        expect(a.indexWhere((e) => e == 99), -1);
        expect(a.indexWhere((e) => e == 2, 2), 2);
        expect(a.indexWhere((e) => e == 2, 3), -1);
      });

      test('insert', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(() => a.insert(1, 5), throwsUnsupportedError);
      });

      test('insertAll', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(() => a.insertAll(1, [5, 6, 7, 8]), throwsUnsupportedError);
      });

      test('emptyness', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.isEmpty, false);
        expect(a.isNotEmpty, true);

        final b = Array.fromList(<int>[]);
        expect(b.isEmpty, true);
        expect(b.isNotEmpty, false);
      });

      test('iterator', () {
        final a = Array.fromList([1, 2, 3, 4]);
        final iter = a.iterator;
        expect([for (; iter.moveNext();) iter.current], [1, 2, 3, 4]);
      });

      test('join', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.join(','), '1,2,3,4');
      });

      test('lastIndexOf', () {
        final a = Array.fromList([1, 2, 3, 4, 1, 2, 3, 4]);
        expect(a.lastIndexOf(3), 7);
        expect(a.lastIndexOf(99), -1);
        expect(a.lastIndexOf(3, 7), 7);
        expect(a.lastIndexOf(3, 6), 3);
        expect(a.lastIndexOf(3, 3), 3);
        expect(a.lastIndexOf(3, 2), -1);
      });

      test('lastIndexWhere', () {
        final a = Array.fromList([1, 2, 3, 4, 1, 2, 3, 4]);
        expect(a.lastIndexWhere((e) => e == 3), 7);
        expect(a.lastIndexWhere((e) => e == 99), -1);
        expect(a.lastIndexWhere((e) => e == 3, 7), 7);
        expect(a.lastIndexWhere((e) => e == 3, 6), 3);
        expect(a.lastIndexWhere((e) => e == 3, 3), 3);
        expect(a.lastIndexWhere((e) => e == 3, 2), -1);
      });

      test('lastWhere', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.lastWhere((e) => e % 2 == 0), 4);
        expect(a.lastWhere((e) => e % 5 == 0, orElse: () => 99), 99);
      });

      test('map', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.map((e) => e % 2), [1, 0, 1, 0]);
      });

      test('reduce', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.reduce((a, b) => a + b), 10);
      });

      test('remove', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(() => a.remove(3), throwsUnsupportedError);
      });

      test('removeAt', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(() => a.removeAt(3), throwsUnsupportedError);
      });

      test('removeLast', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(() => a.removeLast(), throwsUnsupportedError);
      });

      test('removeRange', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(() => a.removeRange(2, 4), throwsUnsupportedError);
      });

      test('removeWhere', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(() => a.removeWhere((e) => e % 2 == 0), throwsUnsupportedError);
      });

      test('replaceRange', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(() => a.replaceRange(2, 4, [99]), throwsUnsupportedError);
      });

      test('retainWhere', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(() => a.retainWhere((e) => e % 2 == 0), throwsUnsupportedError);
      });

      test('reversed', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.reversed, [4, 3, 2, 1]);
      });

      test('setAll', () {
        final a = Array.fromList([1, 2, 3, 4]);
        a.setAll(2, [99]);
        expect(a, [1, 99, 3, 4]);

        final b = Array.fromList([1, 2, 3, 4]);
        expect(() => b.setAll(2, [99, 99, 99, 99]), throwsRangeError);
      });

      test('setRange', () {
        final a = Array.fromList([1, 2, 3, 4]);
        a.setRange(2, 4, [88, 99]); // exact
        expect(a, [1, 88, 99, 4]);

        a.setRange(2, 4, [11, 22, 33]); // bigger
        expect(a, [1, 11, 22, 4]);

        expect(() => a.setRange(2, 4, [999]), throwsStateError); // smaller
      });

      test('shuffle', () {
        final a = Array.fromList([1, 2, 3, 4]);
        a.shuffle();
        expect(a, containsAll([1, 2, 3, 4]));
      });

      test('single', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(() => a.single, throwsStateError);

        final b = Array.fromList(<int>[]);
        expect(() => b.single, throwsStateError);

        final c = Array.fromList([999]);
        expect(c.single, 999);
      });

      test('singleWhere', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(() => a.singleWhere((e) => e % 2 == 1), throwsStateError);

        expect(() => a.singleWhere((e) => e % 5 == 0), throwsStateError);

        expect(a.singleWhere((e) => e % 3 == 0), 3);
      });

      test('skip', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.skip(2), [3, 4]);
      });

      test('skipWhile', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.skipWhile((e) => e == 1), [2, 3, 4]);
      });

      test('sort', () {
        final a = Array.fromList([3, 4, 1, 2]);
        a.sort();
        expect(a, [1, 2, 3, 4]);
      });

      test('sublist', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.sublist(2), [2, 3, 4]);
        final b = a.sublist(2, 4);
        expect(b, [2, 3]);
        b[0] = 999; // check if it's a copy
        expect(a, [1, 2, 3, 4]);
      });

      test('take', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.take(2), [1, 2]);
      });

      test('takeWhile', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.takeWhile((e) => e < 3), [1, 2]);
      });

      test('toList', () {
        final a = Array.fromList([1, 2, 3, 4]);
        final b = a.toList();
        expect(b, [1, 2, 3, 4]);
        b[0] = 999; // check if it's a copy
        expect(a, [1, 2, 3, 4]);
      });

      test('toSet', () {
        final a = Array.fromList([1, 1, 2, 2, 3, 3, 4, 4]);
        expect(a.toSet(), [1, 2, 3, 4]);
      });

      test('where', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.where((e) => e % 2 == 0), [2, 4]);
      });

      test('whereType', () {
        final a = Array.fromList([1, 2, 3, 4]);
        expect(a.whereType<int>(), [1, 2, 3, 4]);
        expect(a.whereType<double>(), <double>[]);
      });
    });
  });

  group('Matrix', () {
    test('fromList', () {
      final m = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ]);
      expect(m.toData(), [1, 4, 7, 2, 5, 8, 3, 6, 9]);
      expect(m[1][1], 1);
      expect(m[1][2], 2);
      expect(m[1][3], 3);
      expect(m[2][1], 4);
      expect(m[2][2], 5);
      expect(m[2][3], 6);
      expect(m[3][1], 7);
      expect(m[3][2], 8);
      expect(m[3][3], 9);
    });

    test('copy', () {
      final m = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ]);
      final n = m.copy();
      expect(n.toData(), [1, 4, 7, 2, 5, 8, 3, 6, 9]);
      n[1][1] = 999;
      expect(n.toData(), [999, 4, 7, 2, 5, 8, 3, 6, 9]);
      expect(m.toData(), [1, 4, 7, 2, 5, 8, 3, 6, 9]);
      final o = m(2, 2);
      expect(o.toData().sublist(0, 4), [5, 8, 3, 6]);
      final p = o.copy();
      expect(p.toData().sublist(0, 4), [5, 8, 3, 6]);
    });

    test('non squared horizontal', () {
      final m = Matrix.fromList([
        [1, 2, 3, 4, 5, 6],
        [7, 8, 9, 10, 11, 12],
      ]);
      expect(m.toData(), [1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12]);
      expect([for (var i = 1; i <= 6; i++) m[1][i]], [1, 2, 3, 4, 5, 6]);
      expect([for (var i = 1; i <= 6; i++) m[2][i]], [7, 8, 9, 10, 11, 12]);
    });

    test('non squared vertical', () {
      final m = Matrix.fromList([
        [1, 7],
        [2, 8],
        [3, 9],
        [4, 10],
        [5, 11],
        [6, 12],
      ]);
      expect(m.toData(), [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]);
      expect([m[1][1], m[1][2]], [1, 7]);
      expect([m[2][1], m[2][2]], [2, 8]);
      expect([m[3][1], m[3][2]], [3, 9]);
      expect([m[4][1], m[4][2]], [4, 10]);
      expect([m[5][1], m[5][2]], [5, 11]);
      expect([m[6][1], m[6][2]], [6, 12]);
    });

    test('offset', () {
      final m = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ], offset: (
        x: -1,
        y: -10
      ));
      expect(m.toData(), [1, 4, 7, 2, 5, 8, 3, 6, 9]);
      expect(m[10][1], 1);
      expect(m[10][2], 2);
      expect(m[10][3], 3);
      expect(m[11][1], 4);
      expect(m[11][2], 5);
      expect(m[11][3], 6);
      expect(m[12][1], 7);
      expect(m[12][2], 8);
      expect(m[12][3], 9);
    });

    test('offset (non squared)', () {
      final m = Matrix.fromList([
        [1, 2, 3, 4, 5, 6],
        [7, 8, 9, 10, 11, 12],
      ], offset: (
        x: -100,
        y: -10
      ));
      expect(m.toData(), [1, 7, 2, 8, 3, 9, 4, 10, 5, 11, 6, 12]);
      expect([for (var i = 100; i < 106; i++) m[10][i]], [1, 2, 3, 4, 5, 6]);
      expect([for (var i = 100; i < 106; i++) m[11][i]], [7, 8, 9, 10, 11, 12]);
    });

    test('assignment', () {
      final m = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ]);
      m[2][3] = 999;
      expect(m[1][1], 1);
      expect(m[1][2], 2);
      expect(m[1][3], 3);
      expect(m[2][1], 4);
      expect(m[2][2], 5);
      expect(m[2][3], 999);
      expect(m[3][1], 7);
      expect(m[3][2], 8);
      expect(m[3][3], 9);
    });

    test('1x1', () {
      final m = Matrix.fromList([
        [1],
      ]);
      expect(m[1][1], 1);
      expect(m.first, 1);

      final a = Array.fromList([1]);
      final m2 = a.asMatrix(30);
      expect(m2[1][1], 1);
      expect(m2.first, 1);
    });

    test('first', () {
      final m = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ]);
      expect(m.first, 1);
      m.first = 999;
      expect(m[1][1], 999);
      expect(m[1][2], 2);
      expect(m[1][3], 3);
      expect(m[2][1], 4);
      expect(m[2][2], 5);
      expect(m[2][3], 6);
      expect(m[3][1], 7);
      expect(m[3][2], 8);
      expect(m[3][3], 9);
    });

    test('sliced first', () {
      final m = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ]);
      final s1 = m(2, 1);
      s1.first = 999;
      expect(s1[1][1], 999);
      expect(s1[1][2], 5);
      expect(s1[1][3], 6);
      expect(s1[2][1], 7);
      expect(s1[2][2], 8);
      expect(s1[2][3], 9);
    });

    test('offset (zero indexed)', () {
      final m = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ], offset: zeroIndexedMatrixOffset);
      expect(m[0][0], 1);
      expect(m[0][1], 2);
      expect(m[0][2], 3);
      expect(m[1][0], 4);
      expect(m[1][1], 5);
      expect(m[1][2], 6);
      expect(m[2][0], 7);
      expect(m[2][1], 8);
      expect(m[2][2], 9);
    });

    test('slice (zero indexed)', () {
      final m = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ], offset: zeroIndexedMatrixOffset);

      final n = m(1, 1);
      expect(n[0][0], 5);
      expect(n[0][1], 6);
      expect(n[1][0], 8);
      expect(n[1][1], 9);

      // back to one-indexed
      final o = m(1, 1, offset: oneIndexedMatrixOffset);
      expect(o[1][1], 5);
      expect(o[1][2], 6);
      expect(o[2][1], 8);
      expect(o[2][2], 9);
    });

    test('sliced assignment', () {
      final m = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ]);

      final s = m(2, 1);
      s[2][3] = 999;
      expect(s[1][1], 4);
      expect(s[1][2], 5);
      expect(s[1][3], 6);
      expect(s[2][1], 7);
      expect(s[2][2], 8);
      expect(s[2][3], 999);
    });

    test('slice line', () {
      final m = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ])(2, 1);
      expect(m[1][1], 4);
      expect(m[1][2], 5);
      expect(m[1][3], 6);
      expect(m[2][1], 7);
      expect(m[2][2], 8);
      expect(m[2][3], 9);
    });

    test('slice column', () {
      final m = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ])(1, 2);
      expect(m[1][1], 2);
      expect(m[1][2], 3);
      expect(m[2][1], 5);
      expect(m[2][2], 6);
      expect(m[3][1], 8);
      expect(m[3][2], 9);
    });

    test('slice LC', () {
      final m = Matrix.fromList([
        [1, 2, 3, 4],
        [5, 6, 7, 8],
        [9, 10, 11, 12],
      ])(2, 2);
      expect(m[1][1], 6);
      expect(m[1][2], 7);
      expect(m[1][3], 8);
      expect(m[2][1], 10);
      expect(m[2][2], 11);
      expect(m[2][3], 12);
    });

    test('box', () {
      final m = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ]);
      expect(m[2][3], 6);
      final box = m.box(2, 3);
      box.value = 999;
      expect(m[1][1], 1);
      expect(m[1][2], 2);
      expect(m[1][3], 3);
      expect(m[2][1], 4);
      expect(m[2][2], 5);
      expect(m[2][3], 999);
      expect(m[3][1], 7);
      expect(m[3][2], 8);
      expect(m[3][3], 9);
    });

    test('sliced box', () {
      final m = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ])(2, 1);
      expect(m[2][3], 9);
      final box = m.box(2, 3);
      box.value = 999;
      expect(m[1][1], 4);
      expect(m[1][2], 5);
      expect(m[1][3], 6);
      expect(m[2][1], 7);
      expect(m[2][2], 8);
      expect(m[2][3], 999);
    });

    test('asArray', () {
      final m = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ]);
      final a = m.asArray();
      expect(a[1], 1);
      expect(a[2], 4);
      expect(a[3], 7);
      expect(a[4], 2);
      expect(a[5], 5);
      expect(a[6], 8);
      expect(a[7], 3);
      expect(a[8], 6);
      expect(a[9], 9);
    });

    test('sliced asArray', () {
      final m = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ]);
      final a = m(2, 1).asArray();
      expect(a[1], 4);
      expect(a[2], 7);
      expect(a[3], 2);
      expect(a[4], 5);
      expect(a[5], 8);
      expect(a[6], 3);
      expect(a[7], 6);
      expect(a[8], 9);
      // expect(a[9], 0);
    });

    test('sliced asArray', () {
      final m = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ]);
      final a = m(1, 2).asArray();
      expect(a[1], 2);
      expect(a[2], 5);
      expect(a[3], 8);
      expect(a[4], 3);
      expect(a[5], 6);
      expect(a[6], 9);
      // expect(a[7], 0);
      // expect(a[8], 0);
      // expect(a[9], 0);
    });

    test('sliced asArray', () {
      final m = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ]);
      final a = m(2, 2).asArray();
      expect(a[1], 5);
      expect(a[2], 8);
      expect(a[3], 3);
      expect(a[4], 6);
      // expect(a[5], 0);
      // expect(a[6], 0);
      // expect(a[7], 0);
      // expect(a[8], 0);
      // expect(a[9], 0);
    });

    test('with new dimension', () {
      final m = Matrix.fromList([
        List.filled(10, 0),
        List.filled(10, 0),
        List.filled(10, 0),
        List.filled(10, 0),
        List.filled(10, 0),
        List.filled(10, 0),
        List.filled(10, 0),
        List.filled(10, 0),
        List.filled(10, 0),
        List.filled(10, 0),
      ]);

      final s = m(2, 7, ld: 3);
      for (var i = 1; i <= 3; i++) {
        for (var j = 1; j <= 3; j++) {
          s[i][j] = i * 10 + j;
        }
      }

      const expected = [
        [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 11, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 21, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 31, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 12, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 22, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 32, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 13, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 23, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 33, 0, 0, 0],
      ];

      final actual = [
        [for (var j = 1; j <= 10; j++) m[1][j]],
        [for (var j = 1; j <= 10; j++) m[2][j]],
        [for (var j = 1; j <= 10; j++) m[3][j]],
        [for (var j = 1; j <= 10; j++) m[4][j]],
        [for (var j = 1; j <= 10; j++) m[5][j]],
        [for (var j = 1; j <= 10; j++) m[6][j]],
        [for (var j = 1; j <= 10; j++) m[7][j]],
        [for (var j = 1; j <= 10; j++) m[8][j]],
        [for (var j = 1; j <= 10; j++) m[9][j]],
        [for (var j = 1; j <= 10; j++) m[10][j]],
      ];

      for (var i = 1; i <= 10; i++) {
        expect(actual[i - 1], expected[i - 1], reason: 'i=$i');
      }
    });

    test('fill sliced', () {
      const NMAX = 132;
      const NEED = 14, //
          M = 10,
          N = 10;
      final m = Matrix<double>(NMAX * NMAX, NEED);
      final s = m.having(ld: NMAX);

      for (var J = 1; J <= N; J++) {
        for (var I = 1; I <= M; I++) {
          s[I][J] = 1;
        }
      }
    });

    test('having', () {
      final m1 = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ]);

      expect(
        m1,
        MatrixEquals([
          [1, 2, 3],
          [4, 5, 6],
          [7, 8, 9],
        ]),
      );

      final m2 = m1.having(offset: zeroIndexedMatrixOffset);
      m2.first = 999;
      expect(m2[0][0], 999);
      expect(m2[0][1], 2);
      expect(m2[1][0], 4);

      expect(m1[1][1], 999); // `having` should get a new reference

      // `ld` less than original
      final m3 = m1.having(ld: 2, offset: zeroIndexedMatrixOffset);

      expect(
        m3,
        MatrixEquals([
          [999, 7, 5, 3],
          [4, 2, 8, 6],
        ]),
      );

      // `ld` greater than original
      final m4 = m1.having(ld: 4);
      expect(m4[1][1], 999);

      expect(
        m4,
        MatrixEquals([
          [999, 5],
          [4, 8],
          [7, 3],
          [2, 6],
        ]),
      );
    });
  });

  group('Matrix3d', () {
    test('fromList', () {
      final m = Matrix3d.fromList([
        [
          [11, 12, 13],
          [14, 15, 16],
          [17, 18, 19],
        ],
        [
          [21, 22, 23],
          [24, 25, 26],
          [27, 28, 29],
        ],
        [
          [31, 32, 33],
          [34, 35, 36],
          [37, 38, 39],
        ],
      ]);
      expect(m.toData(), [
        11, 21, 31, 14, 24, 34, 17, 27, 37, 12, 22, 32, 15, //
        25, 35, 18, 28, 38, 13, 23, 33, 16, 26, 36, 19, 29, 39
      ]);

      for (var i = 1; i <= 3; i++) {
        var b = 1;
        for (var j = 1; j <= 3; j++) {
          for (var k = 1; k <= 3; k++) {
            expect(m[i][j][k], i * 10 + b, reason: 'i=$i, j=$j, k=$k');
            b += 1;
          }
        }
      }
    });

    test('slice', () {
      final m = Matrix3d.fromList([
        [
          [11, 12, 13],
          [14, 15, 16],
          [17, 18, 19],
        ],
        [
          [21, 22, 23],
          [24, 25, 26],
          [27, 28, 29],
        ],
        [
          [31, 32, 33],
          [34, 35, 36],
          [37, 38, 39],
        ],
      ]);
      final s1 = m(3, 2, 1);
      expect([s1[1][1][1], s1[1][1][2], s1[1][1][3]], [34, 35, 36]);
      expect([s1[1][2][1], s1[1][2][2], s1[1][2][3]], [37, 38, 39]);
      final s2 = m(1, 2, 3);
      expect(s2[1][1][1], 16);
      expect(s2[1][2][1], 19);
      expect(s2[2][1][1], 26);
    });

    test('fromList irregular dimesion', () {
      final m = Matrix3d.fromList([
        [
          [11, 12],
          [13, 14],
          [15, 16],
        ],
        [
          [21, 22],
          [23, 24],
          [25, 26],
        ],
        [
          [31, 32],
          [33, 34],
          [35, 36],
        ],
        [
          [41, 42],
          [43, 44],
          [45, 46],
        ],
      ]);

      expect(m.toData(), [
        11, 21, 31, 41, 13, 23, 33, 43, 15, 25, 35, 45, //
        12, 22, 32, 42, 14, 24, 34, 44, 16, 26, 36, 46
      ]);

      expect([m[4][3][1], m[4][3][2]], [45, 46]);
      expect([m[1][1][1], m[1][1][2]], [11, 12]);
      expect([m[2][1][1], m[2][1][2]], [21, 22]);
      expect([m[2][2][1], m[2][2][2]], [23, 24]);
      expect([m[3][2][1], m[3][2][2]], [33, 34]);
      final s = m(3, 2, 1);
      expect([s[1][1][1], s[1][1][2]], [33, 34]);
      expect([s[1][2][1], s[1][2][2]], [35, 36]);
    });
  });
}

class MatrixEquals<T> extends CustomMatcher {
  MatrixEquals(List<List<T>> expected)
      : super('Complex with', 'complex', [
          for (var i = 0; i < expected.length; i++)
            equals(
                [for (var j = 0; j < expected[i].length; j++) expected[i][j]])
        ]);

  @override
  Object? featureValueOf(dynamic actual) => [
        for (var i = 0; i < (actual as Matrix<T>).dimensions.$1; i++)
          [
            for (var j = 0; j < actual.dimensions.$2; j++)
              actual[i - actual.offset.y][j - actual.offset.x]
          ]
      ];
}


/*
The program:
```
subroutine printAsMatrix(D)
  INTEGER D(3, *)
  DO I = 1, 3
    print *, (D(I,J),J=1,3)
  END DO
end

subroutine printAsArray(D)
  INTEGER D(*)
  print *, (D(J),J=1,9)
end

program test
  INTEGER D(9)
  INTEGER E(3, 3)

  D(1)=1
  D(2)=2
  D(3)=3
  D(4)=4
  D(5)=5
  D(6)=6
  D(7)=7
  D(8)=8
  D(9)=9
  PRINT *, (D(J),J=1,9)
  PRINT *, '--------------------'
  call printAsMatrix(D(1))
  PRINT *, '--------------------'
  E(1,1)=1
  E(1,2)=2
  E(1,3)=3
  E(2,1)=4
  E(2,2)=5
  E(2,3)=6
  E(3,1)=7
  E(3,2)=8
  E(3,3)=9
  DO I = 1, 3
    print *, (E(I,J),J=1,3)
  END DO
  PRINT *, '--------------------'
  call printAsArray(E)

end program
```

Outputs:
```
           1           2           3           4           5           6           7           8           9
 --------------------
           1           4           7
           2           5           8
           3           6           9
 --------------------
           1           2           3
           4           5           6
           7           8           9
 --------------------
           1           4           7           2           5           8           3           6           9
```
 */

/*
Sliced
```
subroutine printAsMatrix(D)
  INTEGER D(3, *)
  DO I = 1, 3
    print *, (D(I,J),J=1,3)
  END DO
end

subroutine printAsArray(D)
  INTEGER D(*)
  print *, (D(J),J=1,9)
end

program test
  INTEGER D(9)
  INTEGER E(3, 3)

  D(1)=1
  D(2)=2
  D(3)=3
  D(4)=4
  D(5)=5
  D(6)=6
  D(7)=7
  D(8)=8
  D(9)=9
  PRINT *, (D(J),J=1,9)
  PRINT *, '--------------------'
  call printAsMatrix(D(3))
  PRINT *, '--------------------'
  E(1,1)=1
  E(1,2)=2
  E(1,3)=3
  E(2,1)=4
  E(2,2)=5
  E(2,3)=6
  E(3,1)=7
  E(3,2)=8
  E(3,3)=9
  DO I = 1, 3
    print *, (E(I,J),J=1,3)
  END DO
  PRINT *, '--------------------'
  call printAsArray(E(2,1))

end program
```

Output:
```
           1           2           3           4           5           6           7           8           9
 --------------------
           3           6           9
           4           7           0
           5           8          10
 --------------------
           1           2           3
           4           5           6
           7           8           9
 --------------------
           4           7           2           5           8           3           6           9           0
```
 */

/*
program test
  INTEGER            IPIVOT( 3, 3, 3 )
  DATA               IPIVOT / 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27 /

  DO I = 1, 3
    DO J = 1, 3
      print *, (IPIVOT(I,J, K),K=1,3)
    END DO
  END DO
endprogram

           1          10          19
           4          13          22
           7          16          25
           2          11          20
           5          14          23
           8          17          26
           3          12          21
           6          15          24
           9          18          27
 */
