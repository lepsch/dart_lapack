import 'package:lapack/src/range.dart';
import 'package:test/test.dart';

void main() {
  group('Range', () {
    group('Positive step', () {
      test('half-open range with start < end', () {
        expect(0.to(5), [0, 1, 2, 3, 4]);
      });

      test('closed range with start < end', () {
        expect(0.through(5), [0, 1, 2, 3, 4, 5]);
      });

      test('half-open range with start == end', () {
        expect(5.to(5), <int>[]);
      });

      test('closed range with start == end', () {
        expect(5.through(5), [5]);
      });

      test('half-open range with start > end', () {
        expect(10.to(5), <int>[]);
      });

      test('closed range with start > end', () {
        expect(10.through(5), <int>[]);
      });
    });

    group('Negative step', () {
      test('half-open range with start < end', () {
        expect(0.to(5, step: -1), <int>[]);
      });

      test('closed range with start < end', () {
        expect(0.through(5, step: -1), <int>[]);
      });

      test('half-open range with start == end', () {
        expect(5.to(5, step: -1), <int>[]);
      });

      test('closed range with start == end', () {
        expect(5.through(5, step: -1), [5]);
      });

      test('half-open range with start > end', () {
        expect(10.to(5, step: -1), [10, 9, 8, 7, 6]);
      });

      test('closed range with start > end', () {
        expect(10.through(5, step: -1), [10, 9, 8, 7, 6, 5]);
      });
    });

    group('Step greater than 1', () {
      test('half-open range with start < end', () {
        expect(0.to(5, step: 2), [0, 2, 4]);
      });

      test('closed range with start < end', () {
        expect(0.through(5, step: 2), [0, 2, 4]);
      });

      test('half-open range with start == end', () {
        expect(5.to(5, step: 2), <int>[]);
      });

      test('closed range with start == end', () {
        expect(5.through(5, step: 2), [5]);
      });

      test('half-open range with start > end', () {
        expect(10.to(5, step: 2), <int>[]);
      });

      test('closed range with start > end', () {
        expect(10.through(5, step: 2), <int>[]);
      });
    });

    group('Step lower than -1', () {
      test('half-open range with start < end', () {
        expect(0.to(5, step: -2), <int>[]);
      });

      test('closed range with start < end', () {
        expect(0.through(5, step: -2), <int>[]);
      });

      test('half-open range with start == end', () {
        expect(5.to(5, step: -2), <int>[]);
      });

      test('closed range with start == end', () {
        expect(5.through(5, step: -2), [5]);
      });

      test('half-open range with start > end', () {
        expect(10.to(5, step: -2), [10, 8, 6]);
      });

      test('closed range with start > end', () {
        expect(10.through(5, step: -2), [10, 8, 6]);
      });
    });
  });
}
