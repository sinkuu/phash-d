module phash;


enum int MaxFileSize = (1<<30); /* 1GB file size limit (for mvp files) */
enum size_t HeaderSize = 64;     /* header size for mvp file */


extern (C++) struct CImg(T) { }


extern (C):
nothrow: @nogc:


/* structure for a single hash */
struct DP
{
    char* id;
    void* hash;
    float* path;
	uint hash_length;
    ubyte hash_type;
}

struct BinHash 
{
	ubyte* hash;
	uint bytelength;
	uint byteidx; // used by addbit()
	ubyte bitmask;  // used by addbit()

	/*
	 * add a single bit to hash. the bits are 
	 * written from left to right.
	 */
	int addbit(ubyte bit) pure nothrow @nogc
	{
		if (bitmask == 0) 
		{
			bitmask = 128; // reset bitmask to "10000000"
			byteidx++;     // jump to next byte in array
		}

		if (byteidx >= bytelength) return -1;
		
		if (bit == 1) *(hash + byteidx) |= bitmask;
		bitmask >>=1;
		return 0;
	}	
}

BinHash* _ph_bmb_new(uint bytelength);
void ph_bmb_free(BinHash *binHash);

/*! /brief Radon Projection info
 */
version (pHash_ImageHash)
{
	struct Projections
	{
		CImg!ubyte* R; //contains projections of image of angled lines through center
		int* nb_pix_perline;        //the head of int array denoting the number of pixels of each line
		int size;                   //the size of nb_pix_perline
	}
}

/*! /brief feature vector info
 */
struct Features
{
    double* features;           //the head of the feature array of double's
    int size;                   //the size of the feature array
}


/*! /brief Digest info
 */
struct Digest
{
	char* id;                 //hash id
	ubyte* coeffs;            //the head of the digest integer coefficient array
	int size;                 //the size of the coeff array
}


/* variables for textual hash */
enum int KgramLength = 50;
enum int WindowLength = 100;
enum int delta = 1;

auto ROTATELEFT(X, BITS)(X x, BITS bits)
{
	return (x << bits) | (x >> 64-bits);
}

struct TxtHashPoint
{
	ulong hash;
	size_t index; /*pos of hash in orig file */
}

struct TxtMatch
{
	size_t first_index; /* offset into first file */
	size_t second_index; /* offset into second file */
	uint length;    /*length of match between 2 files */
}

version (pHash_pthread)
{
	int ph_num_threads();
}

/* /brief alloc a single data point
 *  allocates path array, does nto set id or path
 */
DP* ph_malloc_datapoint(int hashtype);

/** /brief free a datapoint and its path
 *
 */
void ph_free_datapoint(DP *dp);

/*! /brief copyright information
 */
const(char*) ph_about() @safe pure;

/*! /brief radon function
 *  Find radon projections of N lines running through the image center for lines angled 0
 *  to 180 degrees from horizontal.
 *  /param img - CImg src image
 *  /param  N  - int number of angled lines to consider.
 *  /param  projs - (out) Projections struct 
 *  /return int value - less than 0 for error
 */
version (pHash_ImageHash)
{
	int ph_radon_projections(ref const(CImg!ubyte) img, int N, ref Projections projs);

	/*! /brief feature vector
	 *         compute the feature vector from a radon projection map.
	 *  /param  projs - Projections struct
	 *  /param  fv    - (out) Features struct
	 *  /return int value - less than 0 for error
	*/
	int ph_feature_vector(ref const Projections projs,ref Features fv);

	/*! /brief dct 
	 *  Compute the dct of a given vector
	 *  /param R - vector of input series
	 *  /param D - (out) the dct of R
	 *  /return  int value - less than 0 for error
	*/
	int ph_dct(ref const Features fv, ref Digest digest);

	/*! /brief cross correlation for 2 series
	 *  Compute the cross correlation of two series vectors
	 *  /param x - Digest struct
	 *  /param y - Digest struct
	 *  /param pcc - double value the peak of cross correlation
	 *  /param threshold - double value for the threshold value for which 2 images
	 *                     are considered the same or different.
	 *  /return - int value - 1 (true) for same, 0 (false) for different, < 0 for error
	 */

	int ph_crosscorr(ref const Digest x,ref const Digest y,ref double pcc, double threshold = 0.90);

	/*! /brief image digest
	 *  Compute the image digest for an image given the input image
	 *  /param img - CImg object representing an input image
	 *  /param sigma - double value for the deviation for a gaussian filter function 
	 *  /param gamma - double value for gamma correction on the input image
	 *  /param digest - (out) Digest struct
	 *  /param N      - int value for the number of angles to consider. 
	 *  /return       - less than 0 for error
	 */
	int _ph_image_digest(ref const CImg!ubyte img,double sigma, double gamma,ref Digest digest,int N=180);

	/*! /brief image digest
	 *  Compute the image digest given the file name.
	 *  /param file - string value for file name of input image.
	 *  /param sigma - double value for the deviation for gaussian filter
	 *  /param gamma - double value for gamma correction on the input image.
	 *  /param digest - Digest struct
	 *  /param N      - int value for number of angles to consider
	 */
	int ph_image_digest(const char *file, double sigma, double gamma, ref Digest digest,int N=180);


	/*! /brief compare 2 images
	 *  /param imA - CImg object of first image 
	 *  /param imB - CImg object of second image
	 *  /param pcc   - (out) double value for peak of cross correlation
	 *  /param sigma - double value for the deviation of gaussian filter
	 *  /param gamma - double value for gamma correction of images
	 *  /param N     - int number for the number of angles of radon projections
	 *  /param theshold - double value for the threshold
	 *  /return int 0 (false) for different images, 1 (true) for same image, less than 0 for error
	 */
	int _ph_compare_images(ref const CImg!ubyte imA,ref const CImg!ubyte imB,ref double pcc, double sigma = 3.5, double gamma = 1.0,int N=180,double threshold=0.90);

	/*! /brief compare 2 images
	 *  Compare 2 images given the file names
	 *  /param file1 - char string of first image file
	 *  /param file2 - char string of second image file
	 *  /param pcc   - (out) double value for peak of cross correlation
	 *  /param sigma - double value for deviation of gaussian filter
	 *  /param gamma - double value for gamma correction of images
	 *  /param N     - int number for number of angles
	 *  /return int 0 (false) for different image, 1 (true) for same images, less than 0 for error
	 */
	int ph_compare_images(const char *file1, const char *file2,ref double pcc, double sigma = 3.5, double gamma=1.0, int N=180,double threshold=0.90);

	/*! /brief return dct matrix, C
	 *  Return DCT matrix of sqare size, N
	 *  /param N - int denoting the size of the square matrix to create.
	 *  /return CImg!double size NxN containing the dct matrix
	 */
	CImg!float* ph_dct_matrix(const int N);

	/*! /brief compute dct robust image hash
	 *  /param file string variable for name of file
	 *  /param hash of type ulong (must be 64-bit variable)
	 *  /return int value - -1 for failure, 1 for success
	 */
	int ph_dct_imagehash(const char* file,ref ulong hash);

	int ph_bmb_imagehash(const char *file, ubyte method, BinHash **ret_hash);
}

version (pHash_pthread)
{
	DP** ph_dct_image_hashes(char*[] files, int count, int threads = 0);
}

version (pHash_VideoHash)
{
	CImgList!ubyte* ph_getKeyFramesFromVideo(const char *filename);

	ulong* ph_dct_videohash(const char *filename, ref int Length);

	DP** ph_dct_video_hashes(char*[] files, int count, int threads = 0);

	double ph_dct_videohash_dist(ulong *hashA, int N1, ulong *hashB, int N2, int threshold=21);
}

/* ! /brief dct video robust hash
 *   Compute video hash based on the dct of normalized video 32x32x64 cube
 *   /param file name of file
 *   /param hash ulong value for hash value
 *   /return int value - less than 0 for error
 */
version (pHash_ImageHash)
{
	int ph_hamming_distance(const ulong hash1,const ulong hash2);

	/** /brief create a list of datapoint's directly from a directory of image files
	 *  /param dirname - path and name of directory containg all image file names
	 *  /param capacity - int value for upper limit on number of hashes
	 *  /param count - number of hashes created (out param)
	 *  /return pointer to a list of DP pointers (NULL for error)
	 */

	DP** ph_read_imagehashes(const char *dirname,int capacity, ref int count);

	/** /brief create MH image hash for filename image
	*   /param filename - string name of image file
	*   /param N - (out) int value for length of image hash returned
	*   /param alpha - int scale factor for marr wavelet (default=2)
	*   /param lvl   - int level of scale factor (default = 1)
	*   /return ubyte array
	**/
	ubyte* ph_mh_imagehash(const char *filename, ref int N, float alpha=2.0f, float lvl = 1.0f);
}
/** /brief count number bits set in given byte
*   /param val - ubyte byte value
*   /return int value for number of bits set
**/
int ph_bitcount8(ubyte val) @safe pure;

/** /brief compute hamming distance between two byte arrays
 *  /param hashA - byte array for first hash
 *  /param lenA - int length of hashA 
 *  /param hashB - byte array for second hash
 *  /param lenB - int length of hashB
 *  /return double value for normalized hamming distance
 **/
double ph_hammingdistance2(ubyte *hashA, int lenA, ubyte *hashB, int lenB);

/** /brief get all the filenames in specified directory
 *  /param dirname - string value for path and filename
 *  /param cap - int value for upper limit to number of files
 *  /param count - int value for number of file names returned
 *  /return array of pointers to string file names (NULL for error)
 **/

char** ph_readfilenames(const char *dirname,ref int count);


/** /brief textual hash for file
 *  /param filename - char* name of file
 *  /param nbpoints - int length of array of return value (out)
 *  /return TxtHashPoint* array of hash points with respective index into file.
 **/
TxtHashPoint* ph_texthash(const char *filename, int *nbpoints);

/** /brief compare 2 text hashes
 *  /param hash1 -TxtHashPoint
 *  /param N1 - int length of hash1
 *  /param hash2 - TxtHashPoint
 *  /param N2 - int length of hash2
 *  /param nbmatches - int number of matches found (out)
 *  /return TxtMatch* - list of all matches
 **/
TxtMatch* ph_compare_text_hashes(TxtHashPoint *hash1, int N1, TxtHashPoint *hash2, int N2, int *nbmatches);

/* random char mapping for textual hash */

extern (D) immutable ulong[256] textkeys = [
	15498727785010036736UL,
	7275080914684608512UL,
	14445630958268841984UL,
	14728618948878663680UL,
	16816925489502355456UL,
	3644179549068984320UL,
	6183768379476672512UL,
	14171334718745739264UL,
	5124038997949022208UL,
	10218941994323935232UL,
	8806421233143906304UL,
	11600620999078313984UL,
	6729085808520724480UL,
	9470575193177980928UL,
	17565538031497117696UL,
	16900815933189128192UL,
	11726811544871239680UL,
	13231792875940872192UL,
	2612106097615437824UL,
	11196599515807219712UL,
	300692472869158912UL,
	4480470094610169856UL,
	2531475774624497664UL,
	14834442768343891968UL,
	2890219059826130944UL,
	7396118625003765760UL,
	2394211153875042304UL,
	2007168123001634816UL,
	18426904923984625664UL,
	4026129272715345920UL,
	9461932602286931968UL,
	15478888635285110784UL,
	11301210195989889024UL,
	5460819486846222336UL,
	11760763510454222848UL,
	9671391611782692864UL,
	9104999035915206656UL,
	17944531898520829952UL,
	5395982256818880512UL,
	14229038033864228864UL,
	9716729819135213568UL,
	14202403489962786816UL,
	7382914959232991232UL,
	16445815627655938048UL,
	5226234609431216128UL,
	6501708925610491904UL,
	14899887495725449216UL,
	16953046154302455808UL,
	1286757727841812480UL,
	17511993593340887040UL,
	9702901604990058496UL,
	1587450200710971392UL,
	3545719622831439872UL,
	12234377379614556160UL,
	16421892977644797952UL,
	6435938682657570816UL,
	1183751930908770304UL,
	369360057810288640UL,
	8443106805659205632UL,
	1163912781183844352UL,
	4395489330525634560UL,
	17905039407946137600UL,
	16642801425058889728UL,
	15696699526515523584UL,
	4919114829672742912UL,
	9956820861803560960UL,
	6921347064588664832UL,
	14024113865587949568UL,
	9454608686614839296UL,
	12317329321407545344UL,
	9806407834332561408UL,
	724594440630435840UL,
	8072988737660780544UL,
	17189322793565552640UL,
	17170410068286373888UL,
	13299223355681931264UL,
	5244287645466492928UL,
	13623553490302271488UL,
	11805525436274835456UL,
	6531045381898240000UL,
	12688803018523541504UL,
	3061682967555342336UL,
	8118495582609211392UL,
	16234522641354981376UL,
	15296060347169898496UL,
	6093644486544457728UL,
	4223717250303000576UL,
	16479812286668603392UL,
	6463004544354746368UL,
	12666824055962206208UL,
	17643725067852447744UL,
	10858493883470315520UL,
	12125119390198792192UL,
	15839782419201785856UL,
	8108449336276287488UL,
	17044234219871535104UL,
	7349859215885729792UL,
	15029796409454886912UL,
	12621604020339867648UL,
	16804467902500569088UL,
	8900381657152880640UL,
	3981267780962877440UL,
	17529062343131004928UL,
	16973370403403595776UL,
	2723846500818878464UL,
	16252728346297761792UL,
	11825849685375975424UL,
	7968134154875305984UL,
	11429537762890481664UL,
	5184631047941259264UL,
	14499179536773545984UL,
	5671596707704471552UL,
	8246314024086536192UL,
	4170931045673205760UL,
	3459375275349901312UL,
	5095630297546883072UL,
	10264575540807598080UL,
	7683092525652901888UL,
	3128698510505934848UL,
	16727580085162344448UL,
	1903172507905556480UL,
	2325679513238765568UL,
	9139329894923108352UL,
	14028291906694283264UL,
	18165461932440551424UL,
	17247779239789330432UL,
	12625782052856266752UL,
	7068577074616729600UL,
	13830831575534665728UL,
	6800641999486582784UL,
	5426300911997681664UL,
	4284469158977994752UL,
	10781909780449460224UL,
	4508619181419134976UL,
	2811095488672038912UL,
	13505756289858273280UL,
	2314603454007345152UL,
	14636945174048014336UL,
	3027146371024027648UL,
	13744141225487761408UL,
	1374832156869656576UL,
	17526325907797573632UL,
	968993859482681344UL,
	9621146180956192768UL,
	3250512879761227776UL,
	4428369143422517248UL,
	14716776478503075840UL,
	13515088420568825856UL,
	12111461669075419136UL,
	17845474997598945280UL,
	11795924440611553280UL,
	14014634185570910208UL,
	1724410437128159232UL,
	2488510261825110016UL,
	9596182018555641856UL,
	1443128295859159040UL,
	1289545427904888832UL,
	3775219997702356992UL,
	8511705379065823232UL,
	15120377003439554560UL,
	10575862005778874368UL,
	13938006291063504896UL,
	958102097297932288UL,
	2911027712518782976UL,
	18446625472482639872UL,
	3769197585969971200UL,
	16416784002377056256UL,
	2314484861370368000UL,
	18406142768607920128UL,
	997186299691532288UL,
	16058626086858129408UL,
	1334230851768025088UL,
	76768133779554304UL,
	17027619946340810752UL,
	10955377032724217856UL,
	3327281022130716672UL,
	3009245016053776384UL,
	7225409437517742080UL,
	16842369442699542528UL,
	15120706693719130112UL,
	6624140361407135744UL,
	10191549809601544192UL,
	10688596805580488704UL,
	8348550798535294976UL,
	12680060080016588800UL,
	1838034750426578944UL,
	9791679102984388608UL,
	13969605507921477632UL,
	5613254748128935936UL,
	18303384482050211840UL,
	10643238446241415168UL,
	16189116753907810304UL,
	13794646699404165120UL,
	11601340543539347456UL,
	653400401306976256UL,
	13794528098177253376UL,
	15370538129509318656UL,
	17070184403684032512UL,
	16109012959547621376UL,
	15329936824407687168UL,
	18067370711965499392UL,
	13720894972696199168UL,
	16664167676175712256UL,
	18144138845745053696UL,
	12301770853917392896UL,
	9172800635190378496UL,
	3024675794166218752UL,
	15311015869971169280UL,
	16398210081298055168UL,
	1420301171746144256UL,
	11984978489980747776UL,
	4575606368995639296UL,
	11611850981347688448UL,
	4226831221851684864UL,
	12924157176120868864UL,
	5845166987654725632UL,
	6064865972278263808UL,
	4269092205395705856UL,
	1368028430456586240UL,
	11678120728997134336UL,
	4125732613736366080UL,
	12011266876698001408UL,
	9420493409195393024UL,
	17920379313140531200UL,
	5165863346527797248UL,
	10073893810502369280UL,
	13268163337608232960UL,
	2089657402327564288UL,
	8697334149066784768UL,
	10930432232036237312UL,
	17419594235325186048UL,
	8317960787322732544UL,
	6204583131022884864UL,
	15637017837791346688UL,
	8015355559358234624UL,
	59609911230726144UL,
	6363074407862108160UL,
	11040031362114387968UL,
	15370625789791830016UL,
	4314540415450611712UL,
	12460332533860532224UL,
	8908860206063026176UL,
	8890146784446251008UL,
	5625439441498669056UL,
	13135691436504645632UL,
	3367559886857568256UL,
	11470606437743329280UL,
	753813335073357824UL,
	7636652092253274112UL,
	12838634868199915520UL,
	12431934064070492160UL,
	11762384705989640192UL,
	6403157671188365312UL,
	3405683408146268160UL,
	11236019945420619776UL,
	11569021017716162560UL];

