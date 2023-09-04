

c_struct = '''
    static int flag;
    static int c_threshold;
    static float notchQueue[ChannelMax][cycleLength + 2];
    static float notchData[ChannelMax][cycleLength];
    static float pbuffer[ChannelMax][cycleLength];
    static int index[ChannelMax];
'''
