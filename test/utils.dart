// Copyright (c) 2024 Guilherme Lepsch. All rights reserved. Use of this
// source code is governed by a BSD-style license that can be found in the
// [LICENSE file](https://github.com/lepsch/dart_lapack/blob/main/LICENSE).

import 'package:path/path.dart' as path;
import 'package:stack_trace/stack_trace.dart' as stacktrace;

String currentFilePath() => path.dirname(stacktrace.Frame.caller(1).uri.toFilePath());
