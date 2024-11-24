# dart_lapack

Pure Dart implementation of LAPACK (Linear Algebra PACKage).

## About

Dart LAPACK is a pure Dart implementation manually converted from the original
[LAPACK in Fortran from Netlib](https://github.com/Reference-LAPACK/lapack).
It's based on the fork of [LAPACK ref 50d68905](https://github.com/Reference-LAPACK/lapack/commits/50d689057af506e256243ff521641454b241a43b)
from January 25 of 2024 (Around LAPACK version 3.12.0).

The package contains:

1. The Dart implementation for LAPACK, and its testing suite;
2. The Dart implementation of the Basic Linear Algebra Subprograms (the Level 1,
   2, and 3 BLAS) needed by Dart LAPACK, and its testing suite;
3. Fortran intrisic subroutines converted to Dart;
4. Base implementations for Matrix, Array, and value type Boxing;

Because [Dart only supports double-precision floating-point numbers](https://dart.dev/guides/language/numbers)
natively, the single-precision real and single-precision complex LAPACK routines
are not included in this package.

## Installation

`dart pub add dart_lapack` or `flutter pub add dart_lapack`

## Usage

```dart
import 'package:dart_lapack/lapack.dart';
```

## TODO

- [x] Basic Linear Algebra Subprograms BLAS;
- [x] Test suite for Basic Linear Algebra Subprograms BLAS;
- [x] Double precision real LAPACK routines;
- [x] Test suite for Double precision real LAPACK routines;
- [ ] Double precision complex LAPACK routines;
- [ ] Test suite for Double precision complex LAPACK routines;
- [ ] Dart doc;

## Testing

Dart LAPACK includes the entire test suite from the original sources converted
from Fortran to Dart.

To run the test suite do:

```shell
dart test
```

## Changelog

Refer to the [Changelog](https://github.com/lepsch/dart_lapack/blob/main/CHANGELOG.md) to get all release notes.

## Contributing

Feature requests and bug reports are welcome on
[GitHub](https://github.com/lepsch/dart_lapack/issues).

## License

BSD 3-Clause Open MPI Variant License, see [LICENSE](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).
