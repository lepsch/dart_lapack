import 'package:lapack/src/install/lsame.dart';
import 'package:lapack/src/format_specifiers_extensions.dart';

void main() {
  int I1, I2;

  // Determine the character set.

  I1 = 'A'.codeUnitAt(0);
  I2 = 'a'.codeUnitAt(0);
  if (I2 - I1 == 32) {
    print(' ASCII String    set');
  } else {
    print(' Non-ASCII String    set, IOFF should be ${I2 - I1}');
  }

  // Test lsame.

  if (!lsame('A', 'A')) _print9999('A', 'A');
  if (!lsame('A', 'a')) _print9999('A', 'a');
  if (!lsame('a', 'A')) _print9999('a', 'A');
  if (!lsame('a', 'a')) _print9999('a', 'a');
  if (lsame('A', 'B')) _print9998('A', 'B');
  if (lsame('A', 'b')) _print9998('A', 'b');
  if (lsame('a', 'B')) _print9998('a', 'B');
  if (lsame('a', 'b')) _print9998('a', 'b');
  if (lsame('O', '/')) _print9998('O', '/');
  if (lsame('/', 'O')) _print9998('/', 'O');
  if (lsame('o', '/')) _print9998('o', '/');
  if (lsame('/', 'o')) _print9998('/', 'o');
  print(' Tests completed');
}

void _print9999(String s1, String s2) {
  print(' *** Error:  LSAME( ${s1.a1}, ${s2.a1}) is .FALSE. ');
}

void _print9998(String s1, String s2) {
  print(' *** Error:  LSAME( ${s1.a1}, ${s2.a1}) is .TRUE. ');
}
