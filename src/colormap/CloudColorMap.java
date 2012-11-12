package colormap;
import java.awt.Color;


public class CloudColorMap implements Mapper<Float, Color>
{
	@Override
	public Color map(Float value)
	{
		//Sky:   .7f .8f .9f
		//Cloud:  DC, 0f,.9f
		
		if (value < .2f)
			return new Color(Color.HSBtoRGB(.7f,  0f, .9f));
		
		else if (value > .8f)
			return new Color(Color.HSBtoRGB(.7f, .8f, .9f));
		
		value = (value - .2f) * (1f / .6f);
		
		float hue        = .7f * value + (1f - value) * .7f;
		float saturation = .8f * value + (1f - value) *  0f;
		float brightness = .9f * value + (1f - value) * .9f;
		
		return new Color(Color.HSBtoRGB(hue, saturation, brightness));
	}
}
