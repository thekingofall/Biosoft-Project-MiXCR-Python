"""MiLib deterministic seeds and Apache Commons Math WELL19937c integer RNG.

Integer arithmetic reproduces Java signed overflow and unsigned right shifts.
WELL19937c is derived from Apache Commons Math (Apache License 2.0); read
hashing and wildcard conversion follow MiLib embedded in MiXCR 3.0.13.
"""
from threading import local, Lock

MASK32 = (1 << 32)-1
MASK64 = (1 << 64)-1
_CODES = {c:i for i,c in enumerate("AGCTNRYSWKMBDHV")}
_WILDCARDS = {"N":"AGCT", "R":"AG", "Y":"CT", "S":"GC", "W":"AT",
              "K":"GT", "M":"AC", "B":"GCT", "D":"AGT", "H":"ACT", "V":"AGC"}


def signed32(value):
    value &= MASK32
    return value if value < 1 << 31 else value - (1 << 32)


def signed64(value):
    value &= MASK64
    return value if value < 1 << 63 else value - (1 << 64)


def jenkin_wang32(key):
    key = (~key + (key << 15)) & MASK32
    key ^= key >> 12
    key = (key + (key << 2)) & MASK32
    key ^= key >> 4
    key = (key * 2057) & MASK32
    return signed32(key ^ (key >> 16))


def jenkin_wang64(key):
    key = ((key ^ MASK64) + (key << 21)) & MASK64
    key ^= key >> 24
    key = (key + (key << 3) + (key << 8)) & MASK64
    key ^= key >> 14
    key = (key + (key << 2) + (key << 4)) & MASK64
    key ^= key >> 28
    return signed64(key + (key << 31))


def java_string_hash(text):
    # Java iterates UTF-16 code units, including surrogate pairs.
    units = text.encode("utf-16-be", errors="surrogatepass")
    value = 0
    for i in range(0, len(units), 2):
        value = (31*value + (units[i] << 8) + units[i+1]) & MASK32
    return signed32(value)


def sequence_hash(sequence):
    # NucleotideSequence inherits the AbstractArraySequence override, which
    # uses Arrays.hashCode(byte[]), rather than Sequence's generic hash.
    value = 1
    for base in sequence.upper():
        value = signed32(value*31 + _CODES[base])
    return signed32(31*jenkin_wang32(1) + value)


def replace_wildcards(read_id, sequence, quality=None):
    """Reproduce UnsafeFactory.fastqParse(..., replaceWildcards=true).

    FASTQ description is irrelevant to this conversion; it uses the read ID.
    IUPAC replacement changes the corresponding Phred+33 quality to '!'.
    """
    seed = signed64(read_id)
    bases = list(sequence.upper().replace(".", "N"))
    qualities = list(quality) if quality is not None else None
    if qualities is not None and len(qualities) != len(bases):
        raise ValueError("Sequence and quality lengths differ")
    for i, base in enumerate(bases):
        if base in _WILDCARDS:
            seed = jenkin_wang64(signed64(seed + i))
            selection = jenkin_wang64(seed)
            selection = signed64(-selection) if selection < 0 else selection
            choices = _WILDCARDS[base]
            # Long.MIN_VALUE cannot be made positive by Java negation.
            # Its remainder is zero for the power-of-two wildcard sizes.
            remainder = abs(selection) % len(choices)
            if selection < 0 and remainder:
                raise ValueError("Java wildcard index overflow for Long.MIN_VALUE")
            bases[i] = choices[remainder]
            if qualities is not None:
                qualities[i] = "!"
        elif base not in "AGCT":
            raise ValueError("Invalid nucleotide: " + base)
    return "".join(bases), "".join(qualities) if qualities is not None else None


def read_seed(read_id, name, sequence, quality=None):
    """Single-read VDJCAligner.hash after FASTQ wildcard replacement.

    Pass the description without its '@' marker, retaining whitespace. Null
    descriptions use the read ID fallback; empty descriptions hash normally.
    """
    sequence, _ = replace_wildcards(read_id, sequence, quality)
    value = signed64(31 + sequence_hash(sequence))
    return signed64(31*value + (java_string_hash(name) if name is not None else jenkin_wang64(read_id)))


def paired_read_seed(read_id, reads):
    """Paired-read hash; reads contains (description, sequence[, quality])."""
    value = 1
    for read in reads:
        name, sequence = read[:2]
        quality = read[2] if len(read) > 2 else None
        sequence, _ = replace_wildcards(read_id, sequence, quality)
        value = signed64(31*value + sequence_hash(sequence))
        value = signed64(31*value + (java_string_hash(name) if name is not None else jenkin_wang64(read_id)))
    return value


class Well19937c:
    """Apache Commons Math long-seed WELL19937c and integer draws."""
    def __init__(self, seed):
        self.reseed(seed)

    def reseed(self, seed):
        seed &= MASK64
        return self.seed_array([seed >> 32, seed & MASK32])

    def seed_array(self, seed):
        if not seed:
            raise ValueError("Seed array must not be empty")
        size = 624
        values = [x & MASK32 for x in seed[:size]]
        for i in range(len(values), size):
            previous = signed32(values[i-len(seed)])
            values.append((1812433253*(previous ^ (previous >> 30)) + i) & MASK32)
        self.values = values
        self.index = 0
        return self

    def next_bits(self, bits):
        if not 1 <= bits <= 32:
            raise ValueError("bits must be from 1 to 32")
        v = self.values
        index = self.index
        rm1, rm2 = (index+623)%624, (index+622)%624
        v0, vm1, vm2, vm3 = v[index], v[(index+70)%624], v[(index+179)%624], v[(index+449)%624]
        z0 = (v[rm1] & 0x80000000) ^ (v[rm2] & 0x7fffffff)
        z1 = (v0 ^ (v0 << 25) ^ vm1 ^ (vm1 >> 27)) & MASK32
        z2 = ((vm2 >> 9) ^ vm3 ^ (vm3 >> 1)) & MASK32
        z3 = z1 ^ z2
        z4 = (z0 ^ z1 ^ (z1 << 9) ^ z2 ^ (z2 << 21) ^ z3 ^ (z3 >> 21)) & MASK32
        v[index], v[rm1], v[rm2] = z3, z4, v[rm2] & 0x80000000
        self.index = rm1
        z4 ^= (z4 << 7) & 0xe46e1700
        z4 ^= (z4 << 15) & 0x9b868000
        return (z4 & MASK32) >> (32-bits)

    def next_int(self, bound=None):
        if bound is None:
            return signed32(self.next_bits(32))
        if type(bound) is not int or not 0 < bound <= 0x7fffffff:
            raise ValueError("bound must be a positive signed Java int")
        if bound & -bound == bound:
            return (bound*self.next_bits(31)) >> 31
        while True:
            bits = self.next_bits(31)
            value = bits % bound
            if signed32(bits-value+bound-1) >= 0:
                return value

    def next_long(self, bound=None):
        if bound is None:
            return signed64((self.next_bits(32) << 32) | self.next_bits(32))
        if type(bound) is not int or not 0 < bound <= 0x7fffffffffffffff:
            raise ValueError("bound must be a positive signed Java long")
        while True:
            bits = (self.next_bits(31) << 32) | self.next_bits(32)
            value = bits % bound
            if signed64(bits-value+bound-1) >= 0:
                return value

    def next_double(self):
        return ((self.next_bits(26) << 26) | self.next_bits(26)) * 2.220446049250313e-16

    def next_boolean(self):
        return bool(self.next_bits(1))

    nextInt = next_int
    nextLong = next_long
    nextDouble = next_double
    nextBoolean = next_boolean
    setSeed = reseed


_thread_local = local()
_seed_counter = 641
_seed_lock = Lock()


def set_global_initial_seed(seed):
    global _seed_counter
    with _seed_lock:
        _seed_counter = signed64(seed)


def get_thread_local_random():
    global _seed_counter
    if not hasattr(_thread_local, "generator"):
        with _seed_lock:
            _seed_counter = signed64(_seed_counter + 353)
            seed = jenkin_wang64(_seed_counter)
        _thread_local.generator = Well19937c(seed)
    return _thread_local.generator


def reseed_thread_local(seed=None):
    generator = get_thread_local_random()
    seed = generator.next_long() if seed is None else seed
    generator.reseed(seed)
    return seed
