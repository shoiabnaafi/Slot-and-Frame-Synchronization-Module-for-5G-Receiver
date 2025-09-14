# Slot-and-Frame-Synchronization-Module-for-5G-Receiver
Main Takeaway: This project implements a low-latency, cross-correlation based slot and frame synchronization module for 5G NR physical layer, demonstrating expertise in timing acquisition, sliding-window algorithms, and real-time C programming for cellular software roles.

Project Overview
This C project locates the start of a 5G NR slot and frame in a continuous IQ sample stream by detecting the Primary Synchronization Signal (PSS) and Secondary Synchronization Signal (SSS). It is organized for clarity, modularity, and real-time constraints, ideal for demonstrating to Apple’s interviewers.

Key Features

PSS correlation detector using optimized FFT-based cross-correlation

SSS correlation detector for frame boundary identification

Adaptive thresholding and peak detection

Timing offset estimation and correction

Sliding-window processing for low latency

Fixed-point and SIMD-friendly code paths
