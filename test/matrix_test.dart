import 'package:lapack/src/complex.dart';
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
      expect(a.toData(), [1, 2, 3, 4]);
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
      expect(a.maxloc(2, 5), 5);
      expect(a.maxloc(2, 4), 4);
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

      final s2 = m[2];
      s2.first = 123;
      expect(s2[1], 123);
      expect(s2[2], 5);
      expect(s2[3], 6);
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
      final s = m.dim(NMAX);

      for (var J = 1; J <= N; J++) {
        for (var I = 1; I <= M; I++) {
          s[I][J] = 1;
        }
      }
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
