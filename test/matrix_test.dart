import 'package:lapack/src/matrix.dart';
import 'package:test/test.dart';

void main() {
  group('Array', () {
    test('fromList', () {
      final a = Array.fromList([1, 2, 3, 4]);
      expect(a[1], 1);
      expect(a[2], 2);
      expect(a[3], 3);
      expect(a[4], 4);
      expect(a.toRawList(), [1, 2, 3, 4]);
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
      expect(a.toRawList(), [2, 3, 4]);

      final b = a(2);
      expect(b[1], 3);
      expect(b[2], 4);
      expect(b.toRawList(), [3, 4]);
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
  });

  group('Matrix', () {
    test('fromList', () {
      final a = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ]);
      expect(a.toRawList(), [1, 4, 7, 2, 5, 8, 3, 6, 9]);
      expect(a[1][1], 1);
      expect(a[1][2], 2);
      expect(a[1][3], 3);
      expect(a[2][1], 4);
      expect(a[2][2], 5);
      expect(a[2][3], 6);
      expect(a[3][1], 7);
      expect(a[3][2], 8);
      expect(a[3][3], 9);
    });

    test('assignment', () {
      final a = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ]);
      a[2][3] = 999;
      expect(a[1][1], 1);
      expect(a[1][2], 2);
      expect(a[1][3], 3);
      expect(a[2][1], 4);
      expect(a[2][2], 5);
      expect(a[2][3], 999);
      expect(a[3][1], 7);
      expect(a[3][2], 8);
      expect(a[3][3], 9);
    });

    test('sliced assignment', () {
      final a = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ])(2, 1);
      a[2][3] = 999;
      expect(a[1][1], 4);
      expect(a[1][2], 5);
      expect(a[1][3], 6);
      expect(a[2][1], 7);
      expect(a[2][2], 8);
      expect(a[2][3], 999);
    });

    test('slice line', () {
      final a = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ])(2, 1);
      expect(a[1][1], 4);
      expect(a[1][2], 5);
      expect(a[1][3], 6);
      expect(a[2][1], 7);
      expect(a[2][2], 8);
      expect(a[2][3], 9);
    });

    test('slice column', () {
      final a = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ])(1, 2);
      expect(a[1][1], 2);
      expect(a[1][2], 3);
      expect(a[2][1], 5);
      expect(a[2][2], 6);
      expect(a[3][1], 8);
      expect(a[3][2], 9);
    });

    test('slice LC', () {
      final a = Matrix.fromList([
        [1, 2, 3, 4],
        [5, 6, 7, 8],
        [9, 10, 11, 12],
      ])(2, 2);
      expect(a[1][1], 6);
      expect(a[1][2], 7);
      expect(a[1][3], 8);
      expect(a[2][1], 10);
      expect(a[2][2], 11);
      expect(a[2][3], 12);
    });

    test('box', () {
      final a = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ]);
      expect(a[2][3], 6);
      final box = a.box(2, 3);
      box.value = 999;
      expect(a[1][1], 1);
      expect(a[1][2], 2);
      expect(a[1][3], 3);
      expect(a[2][1], 4);
      expect(a[2][2], 5);
      expect(a[2][3], 999);
      expect(a[3][1], 7);
      expect(a[3][2], 8);
      expect(a[3][3], 9);
    });

    test('sliced box', () {
      final a = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ])(2, 1);
      expect(a[2][3], 9);
      final box = a.box(2, 3);
      box.value = 999;
      expect(a[1][1], 4);
      expect(a[1][2], 5);
      expect(a[1][3], 6);
      expect(a[2][1], 7);
      expect(a[2][2], 8);
      expect(a[2][3], 999);
    });
  });
}
