import 'package:lapack/src/format_extensions.dart';
import 'package:lapack/src/matrix.dart';

extension DoubleMatrixExtension on Matrix<double> {
  void debug(String name, int m, int n) {
    print('$name(${m.i6}, ${n.i6}) =');
    for (var i = 1; i <= m; i++) {
      var line = '';
      for (var j = 1; j <= m; j++) {
        final double v = this[i][j] == 0 ? 0 : this[i][j];
        line += v.d10_3;
      }
      print(line);
    }
  }
}

extension DoubleArrayExtension on Array<double> {
  void debug(String name, int s) {
    var line = '$name(${s.i6}) = ';
    for (var i = 1; i <= s; i++) {
      final double v = this[i] == 0 ? 0 : this[i];
      line += v.d10_3;
    }
    print(line);
  }
}

extension BoolArrayExtension on Array<bool> {
  void debug(String name, int s) {
    var line = '$name(${s.i6}) = ';
    for (var i = 1; i <= s; i++) {
      line += (this[i] ? 1 : 0).i2;
    }
    print(line);
  }
}
