#include <iostream>
#include <gmpxx.h>
#include <string>

std::string p_str =   "115792089237316195423570985008687907853269984665640564039457584007913111864739";
std::string a_str =   "115792089237316195423570985008687907853269984665640564039457584007913111752419";
std::string b_str =   "13602384";

std::string x_str =   "44328971593885937857970623207174810055095945000614270339392047863929064377300";
std::string y_str =   "73987224069968535275377617159869580030126023743076722472100521420353122284142";

std::string a_th_str = "8";
std::string d_th_str = "48";

struct Point{
    mpz_class x, y, z;

    Point(){
        x = 0;
        y = 0;
        z = 0;
    }

    void set (std::string x_str, std::string y_str, std::string z_str){
        x = x_str;
        y = y_str;
        z = z_str;
    }

    void show(){
        std::cout << "x= " << x << "\n";
        std::cout << "y= " << y << "\n";
        std::cout << "z= " << z << "\n\n";

    }

    void show_aff(){
        std::cout << "x= " << x << "\n";
        std::cout << "y= " << y << "\n\n";
    }


};


int is_equial_aff(Point p1, Point p2){
    if (p1.x == p2.x and p1.y == p2.y){
        return 1;
    }
    return 0;
}

struct params_of_weierstass{
    mpz_class a, b, p, m;
    Point point, result_point;

    params_of_weierstass(){
        a = 0; b = 0; p = 0;
    }

    void set_point(std::string x_str, std::string y_str, std::string z_str){
        point.set(x_str, y_str, z_str);
    }

    params_of_weierstass(std::string a_str, std::string  b_str, std::string p_str){
        a = a_str;
        b = b_str;
        p = p_str;
    }

    void show(){
        std::cout << "PARAMETRS OF WIER:\n";
        std::cout << "a: " << a << "\n";
        std::cout << "b: " << b << "\n";
        std::cout << "p: " << p << "\n";
    }

};


struct twisted_hessian_curve{
    mpz_class a, d, p;
    Point point;

    twisted_hessian_curve(){
        a = 0;
        d = 0;
        p = 0;
    }
    twisted_hessian_curve(std::string a_str, std::string d_str, mpz_class p_){
        a = a_str;
        d = d_str;
        p = p_;
    }


    void point_to_th(struct Point start_point){
        mpz_class buf, buf_x, buf_y,
                invrt, res_x, res_y;
        mpz_class buf_res;

        // formula
        // weierstrass -> twisted hessian curve
        // th_x = (18*d^2 + 72*x)   /  (d^3-12*d*x-108*a+24*y)
        // th_y = (1-(48*y)         /  (d^3-12*d*x-108*a+24*y))
        buf = d*d*d - 12*d*start_point.x - 108*a +24*start_point.y;// одинаковый знаменатель
        buf_x = 18*d*d + 72*start_point.x;
        buf_y = 48*start_point.y;

        buf_res = buf % p;
        if (mpz_sgn (buf_res.get_mpz_t()) == -1){
            buf_res = buf_res + p;
        }
        mpz_invert(invrt.get_mpz_t(), buf_res.get_mpz_t() , p.get_mpz_t());

        res_x = buf_x * invrt;
        res_x = res_x % p;

        res_y = buf_y * invrt;
        res_y = res_y % p;
        res_y = 1 - res_y;
        res_y = res_y % p;
        if (mpz_sgn (res_y.get_mpz_t()) == -1){
            res_y = res_y + p;
        }


        point.x = res_x;
        point.y = res_y;
        point.z = 1;
    }


    Point add_aff(Point point1, Point point2){
        Point R;

        mpz_class x3, y3, down, x3_up, y3_up, down_inv;
        down = a * point1.x * point1.y * point2.x * point2.x - point2.y;
        x3_up = point1.x - point1.y*point1.y * point2.x * point2.y;


        y3_up = point1.y * point2.y*point2.y - a*point1.x*point1.x * point2.x;

        mpz_invert(down_inv.get_mpz_t(), down.get_mpz_t(), p.get_mpz_t());
        x3 = (x3_up * down_inv) % p;
        y3 = (y3_up * down_inv) % p;

        if (mpz_sgn (x3.get_mpz_t()) == -1){
            x3 = x3 + p;
        }
        if (mpz_sgn (y3.get_mpz_t()) == -1){
            y3 = y3 + p;
        }
        R.x = x3;
        R.y = y3;
        return R;
    }

    Point crat_aff(mpz_class k){
        Point R, Q;
        Q.set("0","-1","0");

        R.x = point.x;
        R.y = point.y;

        int number_of_bits;

        number_of_bits = static_cast<int>(mpz_sizeinbase(k.get_mpz_t(), 2));
        for (int i = number_of_bits - 1; i>=0 ;--i){
            if (mpz_tstbit(k.get_mpz_t(), i) == 0){
                R = add_aff(R, Q);
                Q = add_aff(Q, Q);
            }
            else{
                Q = add_aff(Q, R);
                R = add_aff(R, R);
            }
        }
        return Q;
    }

    int check_point_aff(Point point_1){
        mpz_class l_side, r_side;
        l_side = a * (point_1.x * point_1.x * point_1.x) + (point_1.y * point_1.y * point_1.y) + 1;
        r_side = d * point_1.x * point_1.y;
        l_side = l_side % p;
        r_side = r_side % p;

        if (l_side == r_side){

            return 1;
        }
        return 0;
    }

    Point invert_point_aff(Point point1){
        Point res_point;
        mpz_class invert_y;

        mpz_invert(invert_y.get_mpz_t(), point1.y.get_mpz_t(), p.get_mpz_t());

        res_point.x = (point1.x * invert_y) % p;
        res_point.y = invert_y % p;

        return res_point;
    }
};



int main()

{
    //
    // tests for parametric
    //
    params_of_weierstass params (a_str, b_str, p_str);
    twisted_hessian_curve th_curve(a_th_str, d_th_str, params.p);

    params.show();
    params.point.set(x_str, y_str, "1");
    params.point.show_aff();

    std::cout << "TWISTED HESSIAN CURVE\n";
    th_curve.point_to_th(params.point);
    std::cout << "a= " <<th_curve.a << "\n";
    std::cout << "d= " <<th_curve.d << "\n";
    std::cout << "p= " <<th_curve.p << "\n";

    //
    // tests for affine
    //
    std::cout << "tests for affine;\n\n";

    Point Q, q_p, p_q, q_q;
    Q.set("0","-1","0");

    q_p = th_curve.add_aff(Q,th_curve.point);
    p_q = th_curve.add_aff(th_curve.point, Q);
    q_q = th_curve.add_aff(Q, Q);

    std::cout << "(x,y) is on curve? ";
    std::cout << th_curve.check_point_aff(th_curve.point) <<"\n";

    std::cout << "\n\nTEST_1. Neutral element\n";
    std::cout << "q = (0, -1) is neutral element? ";
    if (is_equial_aff(q_p, p_q) and is_equial_aff(q_p, th_curve.point)){ // checking q+p=p+q=p
        std::cout << "true\n";
    }else{
        std::cout << "false\n";
    }


    std::cout << "\n\nTEST_2. small k. [k]P. [k]P is on curve. \n";
    mpz_class k; k = 20;
    Point kp;
    kp = th_curve.crat_aff(k);
    std::cout << "k= " << k ;
    std::cout << "\n[k]P is on curve? ";
    if (th_curve.check_point_aff(kp)){
        std::cout << "true\n";
    }else{
        std::cout << "false\n";
    }


    std::cout << "\n\nTEST_3. small k. [k1 + k2]P = [k1]P + [k2]P. \n";
    mpz_class k1,k2,k3; k1= 3; k2 = 5; k3 =k1+k2;
    Point pk1, pk2, pk3, pk1_2;
    std::cout << "k1=" << k1 <<'\n'  ;
    std::cout << "k2=" << k2 <<'\n'  ;
    std::cout << "k3=" << k3 <<'\n'  ;

    pk1 = th_curve.crat_aff(k1);
    pk2 = th_curve.crat_aff(k2);
    pk1_2 = th_curve.add_aff(pk1, pk2);

    pk3 = th_curve.crat_aff(k3);
    std::cout << "is [k1 + k2]P = [k1]P + [k2]P? ";
    if (is_equial_aff(pk1_2, pk3)){
        std::cout << "true\n";
    }else{
        std::cout << "false\n";
    }


    std::cout << "\n\nTEST_4. m - is order. [m]P = q. \n";
    mpz_class m; m = "115792089237316195423570985008687907853279740477758714817704293727807715164245";
    std::cout << "m=" << m <<'\n' << "result_of_kp:\n" ;
    Point mp_aff;
    mp_aff = th_curve.crat_aff(m);
    mp_aff.y = mp_aff.y - th_curve.p;
    std::cout << "is [m]P = q? ";
    if (is_equial_aff(mp_aff, Q)){
        std::cout << "true\n";
    }else{
        std::cout << "false\n";
    }

    std::cout << "\n\nTEST_5. [m + 1]P = P и [m − 1]P = −P. \n";
    std::cout << "m=" << m <<'\n' << "result_of_kp:\n" ;
    Point m_plus_1_p_aff, m_minus_1_p_aff;
    m_plus_1_p_aff = th_curve.crat_aff(m+1);
    m_minus_1_p_aff = th_curve.crat_aff(m-1);
    std::cout << "is [m+1]P = P? ";

    if (is_equial_aff(m_plus_1_p_aff, th_curve.point)){
        std::cout << "true\n";
    }else{
        std::cout << "false\n";
    }

    std::cout << "is [m-1]P = -P? ";
    if (is_equial_aff(m_minus_1_p_aff, th_curve.invert_point_aff(th_curve.point))){
        std::cout << "true\n";
    }else{
        std::cout << "false\n";
    }

    std::cout << "\n\nTEST_6. big k. [k1 + k2]P = [k1]P + [k2]P. \n";
    mpz_class bk1,bk2,bk3, n;
    Point pbk1, pbk2, pbk3, pbk1_2;
    n = "10000000000000000000000000000000000000";
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    mpz_urandomm(bk1.get_mpz_t(), state, n.get_mpz_t());
    mpz_urandomm(bk2.get_mpz_t(), state, n.get_mpz_t());
    bk3 = bk1 + bk2;

    std::cout << "k1=" << bk1 <<'\n'  ;
    std::cout << "k2=" << bk2 <<'\n'  ;
    std::cout << "k3=" << bk3 <<'\n'  ;

    pbk1 = th_curve.crat_aff(bk1);
    pbk2 = th_curve.crat_aff(bk2);
    pbk1_2 = th_curve.add_aff(pbk1, pbk2);
    pbk3 = th_curve.crat_aff(bk3);

    std::cout << "is [k1 + k2]P = [k1]P + [k2]P? ";
    if (is_equial_aff(pbk1_2, pbk3)){
        std::cout << "true\n";
    }else{
        std::cout << "false\n";
    }


    return 0;
}
