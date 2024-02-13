import 'package:lapack/src/intrinsics/exponent.dart';
import 'package:test/test.dart';

void main() {
  group('Exponent', () {
    test('Known', () {
      expect(exponent(1 / 3), -1);
    });

    test('Literal Positive', () {
      expect(exponent(1.23e+10), 34);
      expect(exponent(123e+10), 41);
      expect(exponent(0.123e+10), 31);
      expect(exponent(0.00123e+10), 24);
      expect(exponent(9.99e+10), 37);
      expect(exponent(1.23e+0), 1);
      expect(exponent(1.23e+1), 4);
    });

    test('Literal Negative', () {
      expect(exponent(1.23e-10), -32);
      expect(exponent(123e-10), -26);
      expect(exponent(0.123e-10), -36);
      expect(exponent(0.00123e-10), -42);
      expect(exponent(9.99e-10), -29);
      expect(exponent(1.23e-0), 1);
      expect(exponent(1.23e-1), -3);
    });
  });
}
