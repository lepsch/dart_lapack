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
      ])(2, 1);
      m[2][3] = 999;
      expect(m[1][1], 4);
      expect(m[1][2], 5);
      expect(m[1][3], 6);
      expect(m[2][1], 7);
      expect(m[2][2], 8);
      expect(m[2][3], 999);
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
