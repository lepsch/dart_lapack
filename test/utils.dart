import 'package:path/path.dart' as path;
import 'package:stack_trace/stack_trace.dart' as stacktrace;

String currentFilePath() => path.dirname(stacktrace.Frame.caller(1).uri.toFilePath());
