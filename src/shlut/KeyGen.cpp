// The things in this file will shift to Offline_Evaluator.cpp and Offline_Evaluator.h
#pragma once

#include <vector>
#include "../utils/types.h"

template<class R>
class SomeType {
    R S;
    R t;
    std::vector<R> sigmas;
    std::vector<R> taus0;
    std::vector<R> taus1;
    R gamma;
    int n;
    int m;

public:
    SomeType() = default;

    std::vector<R> evalAll()
    {
        int num = 1;
        std::vector<R> sLeaves(num);
        std::vector<R> newsLeaves(num * 2);
        std::vector<R> tLeaves(num);
        std::vector<R> newtLeaves(num * 2);

        sLeaves[0] = S;
        tLeaves[0] = t;

        for(int i{0}; i < m ; i++)
        {
            for(int j{0}; j < num; j++)
            {
                R value = prg(sLeaves[j]);
                std::vector<R> values = splitResult(value, n);
                newsLeaves[2 * j] = values[0] ^ (tLeaves[j] * sigmas[i]);
                newsLeaves[2 * j + 1] = values[2] ^ (tLeaves[j] * sigmas[i]);
                newtLeaves[2 * j] = values[1] ^ (tLeaves[j] * taus0[i]);
                newtLeaves[2 * j + 1] = values[3] ^ (tLeaves[j] * taus1[i]);
            }
            num *= 2;
            // Resize vectors
            sLeaves.resize(num);
            tLeaves.resize(num);
            for(int j{0}; j < num; j++)
            {
                sLeaves[j] = newsLeaves[j];
                tLeaves[j] = newtLeaves[j];
            }
        }

        return sLeaves;
    } // Check this function, I think there might be errors while while using 2 sources for coding this up
};

template<class R>
std::vector<R> splitResult(R value, int n) {

    //Assuming value is int but will have to change it slightly for Ring

    Ring part0 = (value >> (n + 2));
    Ring part1 = (value >> (n + 1)) & 1;
    Ring part2 = (value & ((1 << (n + 1)) - 2)) >> 1;
    Ring part3 = value & 1;

    return {part0, part1, part2, part3};
}

template<class R>
std::vector<SomeType<R>> GenerateKey(int n, R alpha, R beta)
{
    R Sa00 = random(n);
    R Sb00 = random(n);
    R ta00 = random(1);
    R tb00 = 1 ^ ta00;
    size_t num = 1;
    std::vector<R> sLeaves(num * 2);
    std::vector<R> newsLeaves(num * 4);
    std::vector<R> tLeaves(num * 2);
    std::vector<R> newtLeaves(num * 4);
    sLeaves[0] = Sa00;
    sLeaves[1] = Sb00;
    tLeaves[0] = ta00;
    tLeaves[1] = tb00;
    std::vector<R> z(2, R(0));
    R sigma;
    std::vector<R> tau(2);
    R alpha_digit;
    std::vector<R> sigmas;
    std::vector<R> taus0;
    std::vector<R> taus1;
    int m = 0;

    while(alpha != 0)
    {
        alpha_digit = alpha % R(2);
        m++;

        for(size_t int i = 0; i < num; i++)
        {
            R value = prg(sLeaves[i]);
            std::vector<R> values = splitResult(value, n);
            newsLeaves[2 * i] = values[0];
            newsLeaves[2 * i + 1] = values[2];
            newtLeaves[2 * i] = values[1];
            newtLeaves[2 * i + 1] = values[3];
            z[0] ^= values[0];
            z[1] ^= values[2];

            size_t j = i + num;
            R value = prg(sLeaves[j]);
            std::vector<R> values = splitResult(value, n);
            newsLeaves[2 * j] = values[0];
            newsLeaves[2 * j + 1] = values[2];
            newtLeaves[2 * j] = values[1];
            newtLeaves[2 * j + 1] = values[3];
            z[0] ^= values[0];
            z[1] ^= values[2];
        }

        sigma = z[1 - alpha_digit];
        tau[0] = (z[0] & 1) ^ alpha_digit ^ R(1);
        tau[1] = (z[1] & 1) ^ alpha_digit;
        sigmas.push_back(sigma);
        taus0.push_back(tau[0]);
        taus1.push_back(tau[1]);

        num *= 2;
        // Resize vectors
        sLeaves.resize(num * 2);
        tLeaves.resize(num * 2);
        for(size_t i = 0; i < num; i++)
        {
            sLeaves[i] = newsLeaves[i] ^ (tLeaves[i / 2] * sigma);
            tLeaves[i] = (newsLeaves[i] & 1) ^ (tLeaves[i / 2] * tau[i & 1]);

            size_t j = i + num;
            sLeaves[j] = newsLeaves[j] ^ (tLeaves[j / 2] * sigma); 
            tLeaves[j] = (newsLeaves[j] & 1) ^ (tLeaves[j / 2] * tau[i & 1]);
        }        
        alpha >> 1;
    }

    R gamma = z[alpha_digit] ^ z[2 + alpha_digit] ^ sigma ^ beta;

    // Make the key class and return
    SomeType<R> keya;
    keya.S = Sa00;
    keya.t = ta00;
    keya.sigmas = sigmas;
    keya.taus0 = taus0;
    keya.taus1 = taus1;
    keya.gamma = gamma;
    keya.n = n;
    keya.m = m;

    SomeType<R> keyb;
    keyb.S = Sb00;
    keyb.t = tb00;
    keyb.sigmas = sigmas;
    keyb.taus0 = taus0;
    keyb.taus1 = taus1;
    keyb.gamma = gamma;
    keyb.n = n;
    keyb.m = m;

    return {keya, keyb};
}