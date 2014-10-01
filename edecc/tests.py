import unittest
import time
import random
import string
#from pysodium import crypto_sign, crypto_scalarmult_curve25519, crypto_scalarmult_curve25519_base
# need to figure out how to use relative imports
from elgamal import ElGamal, PrivateKey
import edwards
import inv
import proj
import mont
import xz
import extended
import ed25519
import curve25519
from eddsa import edwardsPrivateKey
from modular import ModInt

class Test(unittest.TestCase):
    """
    Tests each of the representations + original Edwards for functionality
    and time taken by algebraic operations
    """

    def setUp(self):
        self.ed = ed25519.ed25519()
        self.curve25519 = curve25519.curve25519()
        self.projective = proj.projEdwardsCurve(self.ed)
        self.inverted = inv.invEdwardsCurve(self.ed)
        self.extended = extended.extEdwardsCurve(self.ed)
        #self.mont = xz.xzEdwardsCurve(self.curve25519)
        self.n = 10
        self.msg = "hello"

    def get_params(self, group):
        """ Gets ElGamal and ed25519 public andp private keys"""
        self.e0 = edwardsPrivateKey(group)
        self.e1 = self.e0.public_key()
        self.x0 = PrivateKey(group)
        self.y0 = self.x0.public_key()
        self.x1 = PrivateKey(group)
        self.y1 = self.x1.public_key()

    def test_extended(self):
         self.get_params(self.extended)
         #self.basic(self.extended)
         #self.functionality(self.extended)
         #self.timing(self.extended)
         self.data_plotter(self.extended)

    # def test_projective(self):
    #     self.get_params(self.projective)
    #     self.basic(self.projective)
    #     self.functionality(self.projective)
    #     self.timing(self.projective)
    #     self.data_plotter(self.inverted)

    # def test_inverted(self):
    #     self.basic(self.inverted)
    #     self.inverted_special_points()
    #     self.get_params(self.inverted)
    #     self.functionality(self.inverted)
    #     self.timing(self.inverted)
    #     self.data_plotter(self.inverted)

    # def test_edwards(self):
    #     self.get_params(self.ed)
    #     self.edwards_operations()
    #     self.functionality(self.ed)
    #     self.timing(self.ed)
    #     self.data_plotter(self.ed)

    # def test_mont(self):
    #     self.basic_mont(self.mont)
    #     self.get_params(self.mont)
    #     self.functionality(self.mont)
    #     self.timing(self.mont)
    #     self.data_plotter(self.mont)

    # def test_opt(self):
    #     """
    #     Timing tests for python wrapper of c code from the pysodium package.
    #     """
    #     g = self.ed
    #     n = 100
    #     times1, times2 = 0, 0
    #     print("\nTesting scalar multiplication times: ")
    #     for i in range(n):
    #         s = crypto_scalarmult_curve25519_base(randombytes(crypto_scalarmult_BYTES))
    #         r = crypto_scalarmult_curve25519_base(randombytes(crypto_scalarmult_BYTES))
    #         t0 = time.time()
    #         crypto_scalarmult_curve25519(s,r)
    #         times1 += time.time() - t0
    #         t1 = time.time()
    #         crypto_scalarmult_curve25519_base(randombytes(crypto_scalarmult_BYTES))
    #         times2 += time.time() - t1

    #     print("scalarmult (random): ", times1/n, "mult by base: ", times2/n)

    #     times1, times2, times3 = 0, 0, 0
    #     print("\nTesting signing/verification times: ")
    #     for i in range(n):
    #         msg = ''.join([random.choice(string.ascii_letters + string.digits) \
    #             for n in range(numchars)])
    #         pk, sk = crypto_sign_keypair()

    #         t0 = time.time()
    #         signed = crypto_sign(msg,sk)
    #         times1 += time.time() - t0
    #         changed = signed[:crypto_sign_BYTES]+'0'+signed[crypto_sign_BYTES+1:]
    #         print (crypto_sign_open(signed, pk))
    #         t1 = time.time()
    #         crypto_sign_open(changed, pk)
    #         times2 += time.time() - t1
    #         seed = crypto_generichash(msg, outlen=crypto_sign_SEEDBYTES)
    #         pk2, sk2 = crypto_sign_seed_keypair(seed)
    #         t2 = time.time()
    #         pk, sk = crypto_sign_seed_keypair(seed)
    #         times3 += time.time() - t2
    #         assert pk == pk2 and sk == sk2

    #     print("avgtime: signing: ", times1/n, "verification: ", times2/n,
    #         "w/ hash: ", times3/n)

    #     with open(cwrappernummsg.txt, 'w') as file:
    #         nummessages = 100
    #         for n in range(nummessages):
    #             t0 = time.time()
    #             for n in range(i):
    #                 signed = crypto_sign(self.msg,sk)
    #                 changed = signed[:crypto_sign_BYTES]+'0'+signed[crypto_sign_BYTES+1:]
    #                 crypto_sign_open(changed, pk)
    #             time = time.time() - t0
    #             file.write(n + "\t" + time + "\n")

    #     with open(cwrapperlenmsg.txt, 'w') as file:
    #         numchars = 1000
    #         for n in range(numchars):
    #             signed = crypto_sign(self.msg,sk)
    #             changed = signed[:crypto_sign_BYTES]+'0'+signed[crypto_sign_BYTES+1:]
    #             crypto_sign_open(changed, pk)
    #             time = time.time() - t0
    #         file.write(n + "\t" + time + "\n")

    def data_plotter(self, group):
        """
        Writes outputs of timing tests to text files to be plotted by GNUplot.
        """
        # For each of the coordinate systems, num messages vs. time taken to sign/verify

        times1, times2 = 0, 0
        with open(str(group.name)+"nummsg", 'w') as file:
            nummessages = 100
            for n in range(nummessages):
                t0 = time.time()
                for i in range(n):
                    signature = self.e0.sign(self.msg)
                    self.e1.verify(self.msg, signature)
                time1 = time.time() - t0
                file.write(str(n) + "\t" + str(time1) + "\n")

        #how long it takes for message to be signed/verified vs. length of message
        with open(str(group.name)+"lenmsg", 'w') as file:
            numchars = 1000
            for n in range(numchars):
                t0 = time.time()
                signature = self.e0.sign(self.msg)
                self.e1.verify(self.msg, signature)
                time2 = time.time() - t0
                file.write(str(n) + "\t" + str(time2) + "\n")

    def timing(self, group):
        """
        Huge suite of timing tests, for each operation and signing/verification
        n is number of trials
        """
        g = group
        n = 100
        times1 = 0
        times2 = 0
        times3 = 0
        print("\nTesting addition/multiplication/doubling times: ")
        for i in range(n):
            r1 = g.point().random_element()
            r2 = g.point().random_element()
            k = g.secret().random_secret()
            t0 = time.time()
            r1.add(r1, r2)
            times1 += time.time() - t0
            t1 = time.time()
            r2.multiply(r1, k)
            times2 += time.time() - t1
            t2 = time.time()
            r1.double(r1)
            times3 += time.time() - t2
        print("add: ", times1/n, "mult: ", times2/n, "double: ", times3/n)

        times1, times2 = 0
        numchars = 20
        print("\nTesting encoding/decoding times: ")
        for i in range(n):
            str1 = ''.join([random.choice(string.ascii_letters + string.digits)\
                for n in range(1)])
            data = bytes(str1, encoding='utf-8')
            t0 = time.time()
            p = g.point()
            p.encode(data)
            times1 += time.time() - t0
            t1 = time.time()
            p.decode(p)
            times2 += time.time() - t1
        print("avgtime: encoding: ", times1/n, "decoding: ", times2/n)

        times1, times2 = 0
        print("\nTesting signing/verification times: ")
        for i in range(n):
            msg = ''.join([random.choice(string.ascii_letters + string.digits) \
                    for n in range(1)])
            t0 = time.time()
            signature = self.e0.sign(msg)
            times1 += time.time() - t0
            t1 = time.time()
            self.e1.verify(msg, signature)
            times2 += time.time() - t1
        print("avgtime: signing: ", times1/n, "verification: ", times2/n,
            "total: ", times1/n + times2/n)

        times1 = 0
        print("\nTesting key exchange times: ")
        for i in range(n):
            t0 = time.time()
            self.x0.exchange(self.y1)
            times1 += time.time() - t0
        print("avgtime: ", times1/n)

    def functionality(self, group):
        """
        Tests other ECC functionality.
        """
        self.g = group
        self.encoding()
        self.encryption()
        self.exchange()
        self.ecdsa()

    def encoding(self):
        """
        Test point encoding/decoding with random strings.
        """
        elem = self.g.point()
        msg0 = b"" # 0 character encoding: EDGE CASE
        msg1 = b"a" # 1 character encoding
        msg2 = b"Hello" # multiple character encoding
        msg3 = b"abcdefghijklmnopqrstuvxyzab"# 28 characters + 2 padding
        self.assertEqual(elem.decode(elem.encode(msg0)), msg0)
        self.assertEqual(elem.decode(elem.encode(msg1)), msg1)
        self.assertEqual(elem.decode(elem.encode(msg2)), msg2)
        self.assertEqual(elem.decode(elem.encode(msg3)), msg3)

        for i in range(self.n):
            str1 = ''.join([random.choice(string.ascii_letters + string.digits)\
                for n in range(10)])
            data = bytes(str1, encoding='utf-8')
            self.assertTrue(elem.decode(elem.encode(data)) == data)

    def ecdsa(self):
        """
        Test ecdsa functionality with random strings.
        """
        g = self.g
        for i in range(self.n):
            msg = ''.join([random.choice(string.ascii_letters + string.digits) \
                for n in range(20)])
            signature = self.e0.sign(msg)
            self.e1.verify(msg, signature)

    def exchange(self):
        """
        Tests public key exchange.
        """
        for i in range(self.n):
            self.assertEqual(self.x0.exchange(self.y1), self.y1.exchange(self.x0))
            self.assertEqual(self.x1.exchange(self.y0), self.y0.exchange(self.x1))

    def encryption(self):
        for i in range(self.n):
            str1 = ''.join([random.choice(string.ascii_letters + string.digits)\
                for n in range(10)])
            data = bytes(str1, encoding='utf-8')
            encrypted = self.y0.encrypt(data)
            self.assertTrue(self.x0.decrypt(encrypted))

    def basic_mont(self):
        """
        Montgomery addition is not unified, so must be dealt with separately.
        """
        g = self.g
        self.assertTrue(g.is_element(g.generator()))
        self.assertTrue(g.is_element(g.identity()))
        #self.assertTrue(g.add(g.identity(), g.identity()) == g.identity())
        #self.assertTrue(g.add(g.generator(), g.identity()) == g.generator())
        r1 = g.add(g.generator(), g.generator())
        r2 = g.double(g.generator())
        self.assertTrue(g.to_ep(r1) == g.to_ep(r2))
        g.double(g.identity())
        for i in range(10):
            r1 = g.random_element()
            self.assertTrue(g.g.is_element(g.to_ep(r1)))

    def basic(self, group):
        """
        Set of baseline tests for each coordinate system:
        random element generation, addition, doubling, and multiplication
        """
        self.g = group
        g = self.g
        r1 = g.point()
        r2 = g.point()
        r3 = g.point()
        ed = self.ed.point()
        self.assertTrue(r1.generator()._on_curve())
        self.assertTrue(r1.add(r1.identity(), r1.identity()).equal(r1.identity()))
        self.assertTrue(r2.add(r2.generator(), r2.identity()).equal(r2.generator()))
        r3.add(r1.generator(), r1.generator())
        r2.double(r2.generator())
        self.assertTrue(r2.to_ep(r2).equal(r3.to_ep(r3)))

        self.assertTrue(ed.add(ed.generator(), ed.generator()).equal(r3.to_ep(r3)))
        for i in range(10):
            r1 = g.point().random_element()
            self.assertTrue(r1._on_curve())
            self.assertTrue(r2.to_ep(r1)._on_curve())
        self.add()
        self.double()
        self.multiply()

    def add(self):
        """
        Tests the addition operation by creating random elements and checking
        it against the base class's operations.
        """
        g = self.g
        r0 = g.point()
        ref = self.ed.point()
        for i in range(10):
            r1 = g.point().random_element()
            r2 = g.point().random_element()
            self.assertTrue(r1.to_ep(r1.add(r1, r1.identity())).equal(r1.to_ep(r1)))
            self.assertTrue(r1.to_ep(r1)._on_curve() and r2.to_ep(r2)._on_curve())
            self.assertTrue(r0.to_ep(r0.add(r1, r2)).equal(ref.add(r1.to_ep(r1), r2.to_ep(r2))))
            r1.add(r1.generator(), r1)
            self.assertTrue(r1.to_ep(r0.add(r1, r1)).equal(r1.to_ep(r1.double(r1))))

    def inverted_special_points(self):
        """
        Only used for testing inverted coordinate system.
        """
        g = self.g
        a = edwards.EdwardsPoint(self.ed, g.zero, g.one)
        b = edwards.EdwardsPoint(self.ed, g.zero, g.negone)
        c = edwards.EdwardsPoint(self.ed, g.one, g.zero)
        d = edwards.EdwardsPoint(self.ed, g.negone, g.zero)
        a1 = g.s1
        b1 = g.s2
        c1 = g.s3
        d1 = g.s4

        print(type(a))
        r1 = g.point()
        r2 = g.point()
        ref = self.ed.point()
        ref2 = self.ed.point()

        self.assertTrue(r1.to_ep(a1).equal(a) and r2.from_ep(a).equal(a1))
        self.assertTrue(r1.to_ep(b1).equal(b) and r2.from_ep(b).equal(b1))
        self.assertTrue(r1.to_ep(c1).equal(c) and r2.from_ep(c).equal(c1))
        self.assertTrue(r1.to_ep(d1).equal(d) and r2.from_ep(d).equal(d1))

        self.assertTrue(r1.to_ep(r2.add(a1, a1)).equal(ref.add(a, a)))
        self.assertTrue(r1.to_ep(r2.add(b1, b1)).equal(ref.add(b, b)))
        self.assertTrue(r1.to_ep(r2.add(c1, c1)).equal(ref.add(c, c)))
        self.assertTrue(r1.to_ep(r2.add(d1, d1)).equal(ref.add(d, d)))

        self.assertTrue(r1.double(a1).equal(r2.add(a1, a1)))
        self.assertTrue(r1.double(b1).equal(r2.add(b1, b1)))
        self.assertTrue(r1.double(c1).equal(r2.add(c1, c1)))
        self.assertTrue(r1.double(d1).equal(r2.add(d1, d1)))

    def double(self):
        """
        Tests doubling and unification of addition laws (if applicable)
        """
        g = self.g
        ed = self.ed
        ref = ed.point()
        r0 = g.point()
        for i in range(10):
            r1 = g.point().random_element()
            self.assertTrue(r0.to_ep(r0.double(r1)).equal(ref.double(r1.to_ep(r1))))

    def multiply(self):
        """
        Tests multiplication algorithm.
        """
        g = self.g
        ed = self.ed
        for i in range(10):
            k = g.secret()
            ref = ed.point()
            r0 = g.point()
            prod = g.point().multiply(r0.generator(), k)
            newpt = prod.to_ep(prod)
            self.assertTrue((newpt._on_curve()))
            ref.multiply(ed.point().generator(), k)
            self.assertTrue(newpt.equal(ref))

    def clear_denominators(self):
        """
        Tests if the clear denominator function for the inverted and projective
        coordinate systems works properly.
        """
        pass

    def time_add_other(self):
        """
        For extended coordinate system only: use faster unified addition.
        """
        pass

    def triple(self):
        """
        For inverted system only.
        """
        pass

    def edwards_operations(self):
        """
        As the Edwards class does not inherit from itself, its test cases
        are much simpler.
        """
        g = self.ed
        self.assertTrue(self.x1.element._on_curve())
        self.assertTrue(self.x0.element._on_curve())
        for i in range(self.n):
            elem = g.point()
            k = g.secret()
            newelem = elem.add(elem.random_element(), elem.random_element())
            newelem2 = elem.multiply(elem.random_element(), k)

if __name__ == "__main__":
    unittest.main()

