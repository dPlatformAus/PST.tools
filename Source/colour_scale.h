
struct Pixel_Colour {
    unsigned red, green, blue, alpha;
    Pixel_Colour():red(0),green(0),blue(0),alpha(0){};
};

Pixel_Colour get_colour_zero_blue(const unsigned value, const unsigned maximum){
    Pixel_Colour result; //defines an RGB colour for a number between 0 and 1020
    int scaled_value;
    scaled_value = (int)(1020.0*value/maximum);
    result.green=255;
    if (scaled_value<256) result.green=std::min(255,scaled_value);
    if (scaled_value>765) result.green=std::min(255,(255-(scaled_value-765)));
    result.blue=std::min(255,std::max(0,(510-scaled_value)));
    result.red=std::min(255,std::max(0,(scaled_value-510)));
    return result;
}

Pixel_Colour get_colour(const unsigned value, const unsigned maximum){
    Pixel_Colour result; //defines an RGB colour for a number between 0 and 1275
    int scaled_value;
    scaled_value = (int)(1275.0*value/maximum);
    if (scaled_value<256) result.blue=scaled_value;
    else result.blue=std::min(255,std::max(0,(765-scaled_value)));
    result.green=255;
    if (scaled_value<510) result.green=std::min(255,std::max(0,(scaled_value-255)));
    if (scaled_value>1020) result.green=std::min(255,(255-(scaled_value-1020)));
    result.red=std::min(255,std::max(0,(scaled_value-765)));
    return result;
}

void set_pixel(RGBApixel *the_pixel, Pixel_Colour &the_colour){
    the_pixel->Blue = the_colour.blue;
    the_pixel->Green = the_colour.green;
    the_pixel->Red = the_colour.red;
    the_pixel->Alpha = the_colour.alpha;
}

