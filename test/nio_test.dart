import 'package:lapack/src/matrix.dart';
import 'package:lapack/src/nio.dart';
import 'package:test/test.dart';

void main() {
  group('A group of tests', () {
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
    });
  });
}
