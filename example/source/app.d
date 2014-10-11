import std.conv;
import std.exception;
import std.stdio;
import std.string;

import phash;

void main(string[] args)
{
	// print pHash version
	writeln(ph_about().to!string);

	enforce(args.length == 2);

	Digest digest;
	ph_image_digest(args[1].toStringz, 1.0, 1.0, digest);

	writeln("digest size = ", digest.size);
	writeln(digest.coeffs[0 .. digest.size]);
}
