// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:dart_lapack/src/matrix.dart';
import 'package:dart_lapack/src/nio.dart';
import 'package:test/test.dart';

void main() {
  group('NIO', () {
    test('List-Directed I/O', () async {
      final nin = Nin(Stream.fromIterable([
        '''
'quoted1'                  Quoted string
"quoted 2"                 Quoted string with spaces
"quoted '3'"               Mixed quotes
' quoted "4" '             Mixed quotes oposite
'escaped ''5'' quotes'     Escaped quote
F                          Logical false
T                          Logical true
123                        Integer
12.34                      Double
0 1 2 3 5 9                Integer array
0.1 1.2 0.234e2 0.1e-1 5 9 Double array
string T                   Mixed types
ASDF 123
ASDF

'''
            .codeUnits
      ]));

      expect(await nin.readString(), 'quoted1');
      expect(await nin.readString(), 'quoted 2');
      expect(await nin.readString(), "quoted '3'");
      expect(await nin.readString(), ' quoted "4" ');
      expect(await nin.readString(), "escaped '5' quotes");
      expect(await nin.readBool(), false);
      expect(await nin.readBool(), true);
      expect(await nin.readInt(), 123);
      expect(await nin.readDouble(), 12.34);

      final intArray = Array<int>(6);
      await nin.readArray(intArray, 6);
      expect(intArray.toData(), [0, 1, 2, 3, 5, 9]);

      final doubleArray = Array<double>(6);
      await nin.readArray(doubleArray, 6);
      expect(doubleArray.toData(), [0.1, 1.2, 23.4, 0.01, 5, 9]);

      expect(await nin.read2<String, bool>(), ('string', true));
      expect(await nin.read2<String, int?>(), ('ASDF', 123));
      expect(await nin.read2<String, int?>(), ('ASDF', null));
      expect(await nin.read2<String?, int?>(), (null, null));
    });

    test('readMatrix', () async {
      const delta = 0.00001;
      final nin = Nin(Stream.fromIterable([
        '''
   7
  7.81800D-01  5.65700D-01  7.62100D-01  7.43600D-01  2.55300D-01  4.10000D-01
  1.34000D-02
  6.45800D-01  2.66600D-01  5.51000D-01  8.31800D-01  9.27100D-01  6.20900D-01
  7.83900D-01
  1.31600D-01  4.91400D-01  1.77100D-01  1.96400D-01  1.08500D-01  9.27000D-01
  2.24700D-01
  6.41000D-01  4.68900D-01  9.65900D-01  8.88400D-01  3.76900D-01  9.67300D-01
  6.18300D-01
  8.38200D-01  8.74300D-01  4.50700D-01  9.44200D-01  7.75500D-01  9.67600D-01
  7.83100D-01
  3.25900D-01  7.38900D-01  8.30200D-01  4.52100D-01  3.01500D-01  2.13300D-01
  8.43400D-01
  5.24400D-01  5.01600D-01  7.52900D-01  3.83800D-01  8.47900D-01  9.12800D-01
  5.77000D-01
'''
            .codeUnits
      ]));

      final n = await nin.readInt();
      final a = Matrix<double>(n, n);
      await nin.readMatrix(a, n, n);
      expect(a[1][1], closeTo(7.818e-1, delta));
      expect(a[2][1], closeTo(6.458e-1, delta));
      expect(a[7][7], closeTo(5.77e-1, delta));
    });
  });
}
