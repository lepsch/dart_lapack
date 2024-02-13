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

    test('offset (zero indexed)', () {
      final a = Array.fromList([1, 2, 3, 4], offset: 1);
      expect(a[0], 1);
      expect(a[1], 2);
      expect(a[2], 3);
      expect(a[3], 4);
    });
  });

  group('Matrix', () {
    test('fromList', () {
      final m = Matrix.fromList([
        [1, 2, 3],
        [4, 5, 6],
        [7, 8, 9],
      ]);
      expect(m.toRawList(), [1, 4, 7, 2, 5, 8, 3, 6, 9]);
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

    test('with new dimention', () {
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

      for (var i = 1; i <= 10; i++) {
        expect([for (var j = 1; j <= 10; j++) m[i][j]], expected[i - 1],
            reason: 'i=$i');
      }
    });

    test('fill sliced', () {
      const NMAX = 132;
      const NEED = 14, M = 10, N = 10;
      final m = Matrix<double>(NMAX * NMAX, NEED);
      final s = m(1, 7, ld: 132);

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
      expect(m.toRawList(), [
        11, 14, 17, 21, 24, 27, 31, 34, 37, 12, 15, 18, 22, //
        25, 28, 32, 35, 38, 13, 16, 19, 23, 26, 29, 33, 36, 39
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
      final s = m(3, 2, 1);
      expect([s[1][1][1], s[1][1][2], s[1][1][3]], [34, 35, 36]);
      expect([s[1][2][1], s[1][2][2], s[1][2][3]], [37, 38, 39]);
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
